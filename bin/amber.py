'''
Merge ligand and receptor files to write parameter/topology/coordinate files
for Amber simulations.

Created on May 10, 2018

@author: lvotapka
'''

import os, string
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from copy import deepcopy
import parmed
import pytraj
from adv_template import Adv_template


verbose = True

def make_box_info(box_vectors):
  '''Converts a 3x3 matrix of box vectors to x,y,z,alpha,beta,gamma format for 
  CPPTRAJ input.
  Input:
   - box_vectors: A 3x3 matrix representing the triclinic box vectors
  Output:
   - box_info_string: a string formatted for CPPTRAJ 'box' command
  '''
  box_length = box_vectors[0][0].value_in_unit(angstroms) # extract the number without units
  box_info_string = 'x %f y %f z %f alpha 109.4712190 beta 109.4712190 gamma 109.4712190' % (box_length, box_length, box_length)
  return box_info_string

def autoimage_traj(parm_name, trajin_name, trajout_name, box_info, cpptraj_script_location, cpptraj_exe='cpptraj', writing_frames=()):
  '''Runs the CPPTRAJ autoimage command for a triclinic box simulation.
  Input:
   - parm_name: a string representing the filename of a .prmtop or .parm7 AMBER
       parameter/topology file
   - trajin_name: a non-imaged trajectory to load for imaging
   - box_info: a string representing the triclinic box in x,y,z,alpha,beta,gamma
       format
   - cpptraj_script_location: string for the location to write the cpptraj script
   - cpptraj_exe: an optional string representing the OS command to run CPPTRAJ
  Output:
   - None
   '''
  cpptraj_template = '''parm $PARMFILE
box $BOX_INFO
trajin $TRAJIN
autoimage
trajout $TRAJOUT $FRAME_STR
go
quit
'''

  assert len(writing_frames) < 4, 'When writing the autoimaged trajectory, the format of the writing_frames variable must be: (start, stop, offset)'
  if len(writing_frames) == 0: # then write all frames
    frame_str = ''
  elif len(writing_frames) == 1:
    frame_str = 'start %d' % writing_frames[0] # only include the start
  elif len(writing_frames) == 2:
    frame_str = 'start %d stop %d' % (writing_frames[0], writing_frames[1])
  else: # the length is 3
    frame_str = 'start %d stop %d offset %d' % (writing_frames[0], writing_frames[1], writing_frames[2])
  cpptraj_dict = {'PARMFILE':parm_name, 'TRAJIN':trajin_name, 'TRAJOUT':trajout_name, 'BOX_INFO':box_info, 'FRAME_STR':frame_str} # define template dictionary
  cpptraj_script = Adv_template(cpptraj_template, cpptraj_dict) # fill in the values into the template from the dictionary
  extract_file = open(cpptraj_script_location, 'w') # open the script for writing
  extract_file.write(cpptraj_script.get_output()) # write a cpptraj script
  extract_file.close()
  os.system("%s < %s" % (cpptraj_exe, cpptraj_script_location)) # run cpptraj
  return

class AmberSettings():
  def __init__(self):
    self.leap_template = ''
    self.leap_program = 'tleap'

def amber_building(seekrcalc, milestone, amber_settings):
  '''the pre-minimization procedure for amber ff simulations
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to prepare the AMBER system for
   - amber_settings: the AmberSettings() object that defines settings for the
       tleap calculation.
  Output:
   - None
  '''
  i = milestone.index
  building = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'building')
  working_pdb_base = 'holo' #settings['working_pdb_base']
  pdbfile = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'holo_wet.pdb')
  leap_program = amber_settings.leap_program
  if not os.path.exists(pdbfile):
    print "Skipping milestone %d: no holo structure detected at %s." % (i, pdbfile)
    return
  prmtop = os.path.join(building,working_pdb_base+'.parm7')
  inpcrd = os.path.join(building,working_pdb_base+'.rst7')
  newpdb = os.path.join(building,working_pdb_base+'_leap.pdb')
  
  leap_params = { 'HOLO_WET_PDB':pdbfile,
                  'PRMTOP':prmtop,
                  'INPCRD':inpcrd,
                  'LEAP_OUTPUT_PDB':newpdb,
    }
  leap_template = string.Template(amber_settings.leap_template)
  leap_string = leap_template.substitute(leap_params)

  leapfilename = os.path.join(building,'anchor.leap')
  
  if os.path.exists(prmtop) and os.path.exists(inpcrd):
    if verbose: print "Amber Parm files already exist. Skipping build phase to save time."
    milestone.openmm.prmtop_filename = prmtop
    milestone.openmm.inpcrd_filename = inpcrd
    return # save us a little time by only running this when it matters

  leapfile = open(leapfilename,'w')
  leapfile.write(leap_string)
  leapfile.close()
  
  leapcmd = leap_program+' -f '+leapfilename+' > '+os.path.join(building,'leap.out')
  if verbose: print 'running leap using following command:', leapcmd
  
  errcode = os.system(leapcmd)
  # check to make sure everything is ok
  if errcode != 0: # then we hit a problem
    errormsg = "LEaP did not run properly. See %s/leap.out for details" % building
    raise Exception, errormsg
  if (not os.path.exists(prmtop)) or (not os.path.exists(inpcrd)):
    errormsg = "LEaP did not generated expected prmtop & inpcrd files. See %s/leap.out for details" % building
    raise Exception, errormsg
  milestone.openmm.prmtop_filename = prmtop
  milestone.openmm.inpcrd_filename = inpcrd
  return 

def create_restraints(system, topology, constrained_list):
  '''Freeze receptor and ligand molecules, letting water and ions relax.
  Input:
   - system: the OpenMM system object for this milestone
   - topology: the OpenMM topology object of the MD simulation stages
   - constrained_list: a list of integers representing the indices of the atoms
       to constrain.
  Output:
   - None
  '''
  for i, atom in enumerate(topology.atoms()):
    if i in constrained_list:
      system.setParticleMass(i, 0*dalton)
  return
  

def create_simulation(seekrcalc, milestone, box_vectors=None):
  '''create the OpenMM simulation object.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to prepare the AMBER system for
  Output:
   - None
  '''
  if not milestone.openmm.prmtop_filename or not milestone.openmm.inpcrd_filename:
    print "Amber parameter/topology/coordinate files not found for milestone: %d. Skipping." % milestone.index
    return
  prmtop_filename = milestone.openmm.prmtop_filename
  inpcrd_filename = milestone.openmm.inpcrd_filename
  prmtop = AmberPrmtopFile(prmtop_filename)
  inpcrd = AmberInpcrdFile(inpcrd_filename)
  
  system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer) #, constraints=HBonds) # This is fine because h-bonds are always constrained in water!
  integrator = deepcopy(seekrcalc.min_equil.temp_equil_integrator) #LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
  platform = seekrcalc.openmm.platform #Platform.getPlatformByName('CUDA')
  properties = seekrcalc.openmm.properties #{'CudaDeviceIndex':'0', 'CudaPrecision':'double'}
  
  # add restraints
  create_restraints(system, prmtop.topology, seekrcalc.min_equil.constrained) # TODO: clean this up
  barostat = MonteCarloBarostat(1.0*bar, seekrcalc.master_temperature*kelvin, 25)
  system.addForce(barostat)
  
  simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
  simulation.context.setPositions(inpcrd.positions)
  if box_vectors is not None:
    simulation.context.setPeriodicBoxVectors(box_vectors)
  elif inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
  milestone.openmm.system = system
  milestone.openmm.simulation = simulation
  before_min_structure_name = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'min', 'before_min.pdb')
  save_restart(seekrcalc, milestone, pdb_save_filename=before_min_structure_name)
  return
  
# TODO: revisit this function:
# TODO: name seems inappropriate
# TODO: equilibrated.pdb line is hard-coded
def save_restart(seekrcalc, milestone, pdb_save_filename=None, state_file_name=''):
  '''Save an AMBER inpcrd file for easy restart.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to prepare the AMBER system for
  Output:
   - None
  '''
  if not milestone.openmm.prmtop_filename or not milestone.openmm.inpcrd_filename:
    print "Amber parameter/topology/coordinate files not found for milestone: %d. Skipping." % milestone.index
    return
  state = milestone.openmm.simulation.context.getState(getPositions = True, enforcePeriodicBox = False)
  if state_file_name:
    milestone.openmm.simulation.saveState(state_file_name)
  
  positions = state.getPositions()
  box_vectors = state.getPeriodicBoxVectors()
  '''
  amber_parm = parmed.amber.AmberParm(milestone.openmm.prmtop_filename, milestone.openmm.inpcrd_filename)
  amber_parm.positions = positions
  amber_parm.box_vectors = box_vectors
  if not pdb_save_filename:
    rst7_save_filename_not_imaged = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'temp_equil', 'equilibrated_not_imaged.pdb')
    pdb_save_filename_no_epw = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'temp_equil', 'equilibrated_no_epw.pdb')
    pdb_save_filename = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'temp_equil', 'equilibrated.pdb')
    amber_parm.save(rst7_save_filename_not_imaged, overwrite=True)
    if verbose: print "Autoimaging results of temperature equilibration."
    
    cpptraj_script_location = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'temp_equil', 'image_temp_equil.cpptraj')
    box_info = make_box_info(box_vectors)
    autoimage_traj(milestone.openmm.prmtop_filename, rst7_save_filename_not_imaged, pdb_save_filename_no_epw, box_info, cpptraj_script_location=cpptraj_script_location)
    
    equilibrated = mdtraj.load(pdb_save_filename_no_epw, top=milestone.openmm.prmtop_filename)
    amber_parm.positions = equilibrated.xyz
    amber_parm.save(pdb_save_filename)
  '''
  amber_parm = parmed.amber.AmberParm(milestone.openmm.prmtop_filename, milestone.openmm.inpcrd_filename)
  amber_parm.positions = positions
  amber_parm.box = None
  if not pdb_save_filename:
    rst7_save_filename_not_imaged = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'temp_equil', 'equilibrated_not_imaged.rst7')
    pdb_save_filename_not_imaged = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'temp_equil', 'equilibrated_not_imaged.pdb')
    pdb_save_filename = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'temp_equil', 'equilibrated.pdb')
    amber_parm.save(pdb_save_filename_not_imaged, overwrite=True)
    amber_parm.save(rst7_save_filename_not_imaged, overwrite=True)
    if verbose: print "Autoimaging results of temperature equilibration."
    struct = pytraj.load(pdb_save_filename_not_imaged, milestone.openmm.prmtop_filename)
    struct.superpose(ref=0, mask='@CA')
    struct.autoimage()
    pytraj.save(pdb_save_filename, struct, overwrite=True)
    
  else:
    amber_parm.save(pdb_save_filename, overwrite=True)
  return pdb_save_filename
  