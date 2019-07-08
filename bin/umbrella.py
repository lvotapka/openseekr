'''
Created on June 7, 2018

@author: lvotapka
@author: astokely

Functions and objects for running Umbrella sampling
'''

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import os, time, glob, re
from sys import stdout
from seekr import amber
import mdtraj

verbose = True

def load_last_mdtraj_frame(dcd_filename, prmtop_filename, atom_indices=None):
  '''
  This function returns the last frame of a DCD file as an MDtraj object.
  Input:
   - dcd_filename: string that represents the path to the dcd file to load
   - prmtop_filename: string that represents the path to the parm7 file to load
  Output:
   - lastframe: an mdtraj Trajectory object that contains a single frame: the
     last one in the dcd.
  '''
  if atom_indices is not None:
    mytraj_iter = mdtraj.iterload(dcd_filename, top=prmtop_filename, atom_indices=atom_indices)
  else:
    mytraj_iter = mdtraj.iterload(dcd_filename, top=prmtop_filename) # Trajectory object
  
  for frame in mytraj_iter:
    lastframe = frame[-1]
  return lastframe

def create_spherical_forces(seekrcalc, milestone, system):
  '''
  Add the umbrella force: which maintains the ligand on the surface of the 
  spherical milestone.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to run the simulation for
   - system: the OpenMM system object to add the force to
  Output:
   - None
  '''
  new_force = CustomCentroidBondForce(2, '0.5*k*(distance(g1,g2)-radius)^2')
  k = new_force.addGlobalParameter('k', seekrcalc.umbrella_stage.force_constant)
  r0 = new_force.addGlobalParameter('radius', milestone.radius*angstrom)
  g1 = new_force.addGroup(milestone.atom_selection_1)
  g2 = new_force.addGroup(milestone.atom_selection_2)
  if verbose: print "k:", seekrcalc.umbrella_stage.force_constant, "radius:", milestone.radius*angstrom, "g1:", milestone.atom_selection_1, "g2:", milestone.atom_selection_2
  new_force.addBond([g1, g2], [])
  if verbose: print "new_force.getNumGlobalParameters():", new_force.getNumGlobalParameters()
  if verbose: print "new_force.getNumPerBondParameters():", new_force.getNumPerBondParameters()
  return new_force

def create_planar_z_forces(seekrcalc, milestone, system):
  '''
  Add the umbrella force: which maintains the ligand on the surface of the 
  planar milestone.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to run the simulation for
   - system: the OpenMM system object to add the force to
  Output:
   - None
  '''
  new_force = CustomCentroidBondForce(2, '0.5*k*(z2-z1-offset)^2')
  k = new_force.addGlobalParameter('k', seekrcalc.umbrella_stage.force_constant)
  r0 = new_force.addGlobalParameter('offset', milestone.offset*angstrom)
  g1 = new_force.addGroup(milestone.atom_selection_1)
  g2 = new_force.addGroup(milestone.atom_selection_2)
  if verbose: print("k:", seekrcalc.umbrella_stage.force_constant, "offset:", 
                    milestone.offset*angstrom, "g1:", 
                    milestone.atom_selection_1, "g2:", 
                    milestone.atom_selection_2)
  new_force.addBond([g1, g2], [])
  if verbose: print "new_force.getNumGlobalParameters():", new_force.getNumGlobalParameters()
  if verbose: print "new_force.getNumPerBondParameters():", new_force.getNumPerBondParameters()
  return new_force

def create_ellipsoidal_forces(seekrcalc, milestone, system):
  '''
  Add the umbrella force: which maintains the ligand on the surface of the 
  spherical milestone.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to run the simulation for
   - system: the OpenMM system object to add the force to
  Output:
   - None
  '''
  # TODO: code this...
  '''
  new_force = CustomCentroidBondForce(2, '0.5*k*(distance(g1,g2)-nu)^2')
  k = new_force.addGlobalParameter('k', seekrcalc.umbrella_stage.force_constant)
  r0 = new_force.addGlobalParameter('radius', milestone.radius*angstrom)
  g1 = new_force.addGroup(milestone.atom_selection_1)
  g2 = new_force.addGroup(milestone.atom_selection_2)
  '''
  if verbose: print "k:", seekrcalc.umbrella_stage.force_constant, "radius:", milestone.radius*angstrom, "g1:", milestone.atom_selection_1, "g2:", milestone.atom_selection_2
  new_force.addBond([g1, g2], [])
  if verbose: print "new_force.getNumGlobalParameters():", new_force.getNumGlobalParameters()
  if verbose: print "new_force.getNumPerBondParameters():", new_force.getNumPerBondParameters()
  return new_force

def create_rmsd_forces(seekrcalc, milestone, system):
  '''
  Add the umbrella force: which maintains the ligand on the surface of the 
  spherical milestone.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to run the simulation for
   - system: the OpenMM system object to add the force to
  Output:
   - None
  '''
  assert len(milestone.atom_indices1) == len(milestone.atom_indices2), "milestone must have the same number of indices in both molecules for RMSD."
  n = len(milestone.atom_indices1)
  rmsd_list = []
  for i in range(n):
    rmsd_list.append('distance(g%d, g%d)^2' % (i*2+1,i*2+2))
  rmsd_string = '+'.join(rmsd_list)
  force_string = '0.5*k*(%s-radius^2)' % rmsd_string
  print "RMSD umbrella force using string:", force_string
  new_force = CustomCentroidBondForce(2*n, force_string)
  k = new_force.addGlobalParameter('k', seekrcalc.umbrella_stage.force_constant)
  r0 = new_force.addGlobalParameter('radius', milestone.radius*angstrom)
  groups = []
  for i in range(n):
    groups.append(new_force.addGroup(milestone.atom_indices1[i]))
    groups.append(new_force.addGroup(milestone.atom_indices2[i]))
  new_force.addBond(groups, [])
  if verbose: print "k:", seekrcalc.umbrella_stage.force_constant, "radius:", milestone.radius*angstrom, "groups1:", milestone.atom_indices1, "groups2:", milestone.atom_indices2
  if verbose: print "new_force.getNumGlobalParameters():", new_force.getNumGlobalParameters()
  if verbose: print "new_force.getNumPerBondParameters():", new_force.getNumPerBondParameters()
  return new_force

def prep_umbrella_amber(seekrcalc, milestone):
  '''Prepare a system that will use the AMBER forcefield
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to prep the simulation for
  Output:
   - system: the OpenMM system to return from the prmtop and inpcrd inputs
  '''
  prmtop_filename = milestone.openmm.prmtop_filename
  inpcrd_filename = milestone.openmm.inpcrd_filename
  if verbose: print "opening files:", prmtop_filename, inpcrd_filename, pdb_filename
  prmtop = AmberPrmtopFile(prmtop_filename)
  inpcrd = AmberInpcrdFile(inpcrd_filename)
  system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer,
                               constraints=HBonds)
  return system
  
def prep_umbrella_charmm(seekrcalc, milestone):
  '''Prepare a system that will use the CHARMM forcefield
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to prep the simulation for
  Output:
   - system: the OpenMM system to return from the CHARMM inputs
  '''
  psf = CharmmPsfFile(milestone.openmm.psf_filename)

  # Get the coordinates from the PDB
  pdb = PDBFile(milestone.openmm.umbrella_pdb_filename)

  # Load the parameter set.
  params = CharmmParameterSet(milestone.openmm.rtf_filename, 
                              milestone.openmm.par_filename)

  system = psf.createSystem(params, 
                            nonbondedMethod=1.2*nanometer,
                            nonbondedCutoff=None)
  return system

def launch_umbrella_stage(seekrcalc, milestone, box_vectors=None, traj_name='umbrella1.dcd'):
  '''launch an umbrella sampling job.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to run the simulation for
   - box_vectors: 2d numpy array that represents the box vectors. Typically
       extracted from the inpcrd or rst7 files
   - traj_name: the name of this umbrella trajectory
  Output:
   - ending_box_vectors: An array of three vectors that represents the ending
   state of the periodic box. This may not be the same as was started, but could
   have changed through the course of a constant pressure simulation.
   - umbrella_traj: The absolute path to the file of the umbrella trajectory.
  '''
  if seekrcalc.building.ff.lower() == 'amber':
    system = prep_umbrella_amber(seekrcalc, milestone)
  elif seekrcalc.building.ff.lower() == 'charmm':
    system = prep_umbrella_charmm(seekrcalc, milestone)
  else:
    raise Exception, "Forcefield %s not yet implemented in openseekr" % seekrcalc.building.ff
  
  
  pdb_filename = milestone.openmm.umbrella_pdb_filename
  if os.path.exists(pdb_filename):
    pdb = PDBFile(pdb_filename)
    my_positions = pdb.positions
  else: # This will restart a failed umbrella from the last frame
    basename = os.path.basename(pdb_filename)
    no_ext = os.path.splitext(basename)[0] 
    dcd_filename = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'umbrella', '%s.dcd' % no_ext)
    assert os.path.exists(dcd_filename), "Cannot load DCD or PDB file for umbrella stage, none exist:" + dcd_filename
    print "Restarting failed umbrella stage from last frame of DCD file: ", dcd_filename
    last_fwd_frame = load_last_mdtraj_frame(dcd_filename, milestone.openmm.prmtop_filename)
    my_positions = last_fwd_frame.xyz[0]
    
  integrator = LangevinIntegrator(seekrcalc.master_temperature*kelvin, 1/picosecond, 0.002*picoseconds)
  platform = Platform.getPlatformByName('CUDA')
  
  # TODO: change this back
  properties = seekrcalc.openmm.properties
  
  # add restraints
  if milestone.__class__.__name__ in ['Concentric_Spherical_Milestone']:
    new_force = create_spherical_forces(seekrcalc, milestone, system)
  elif milestone.__class__.__name__ in ['Planar_Z_Milestone']:
    new_force = create_planar_z_forces(seekrcalc, milestone, system)
  else:
    raise Exception, "milestone.type = %s not implemented." % milestone.type
  
  system.addForce(new_force)
  if seekrcalc.umbrella_stage.barostat:
    barostat = MonteCarloBarostat(seekrcalc.umbrella_stage.barostat_pressure, seekrcalc.master_temperature*kelvin, seekrcalc.umbrella_stage.barostat_freq)
    system.addForce(barostat)
  
  simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
  simulation.context.setPositions(my_positions)
  simulation.context.setVelocitiesToTemperature(seekrcalc.master_temperature*kelvin)
  if box_vectors:
    simulation.context.setPeriodicBoxVectors(*box_vectors)
  elif inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
  
  if verbose: print "Running energy minimization on milestone:", milestone.index
  simulation.minimizeEnergy()
  
  umbrella_traj = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'umbrella', traj_name)
  simulation.reporters.append(StateDataReporter(stdout, seekrcalc.umbrella_stage.energy_freq, step=True, potentialEnergy=True, temperature=True, volume=True))
  simulation.reporters.append(DCDReporter(umbrella_traj, seekrcalc.umbrella_stage.traj_freq))
  starttime = time.time()
  
  state_filename = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'umbrella', 'backup.state')
  current_step = 0
  while current_step < seekrcalc.umbrella_stage.steps:
  try:
    simulation.saveState(state_filename)
    print "running %d steps" % seekrcalc.umbrella_stage.traj_freq
    simulation.step(seekrcalc.umbrella_stage.traj_freq)
    current_step = current_step + seekrcalc.umbrella_stage.traj_freq
  except ValueError:
    print "Alert! NaN error detected. Restarting from saved state."
    simulation.loadState(state_filename)
    
  #simulation.step(seekrcalc.umbrella_stage.steps) # old way: simulate steps directly
  print "time:", time.time() - starttime, "s"
  end_state = simulation.context.getState(getPositions=True)
  ending_box_vectors = end_state.getPeriodicBoxVectors()
  milestone.openmm.simulation = simulation
  return ending_box_vectors, umbrella_traj

def generate_umbrella_filenames(seekr_calc, milestone):
  umbrella_file_glob = os.path.join(seekr_calc.project.rootdir, milestone.directory, 'md', 'umbrella', 'umbrella*.dcd')
  existing_umbrella_files = glob.glob(umbrella_file_glob)
  if not existing_umbrella_files: # then the directory is empty, we are starting over
    milestone.openmm.umbrella_pdb_filename = os.path.join(seekr_calc.project.rootdir, milestone.directory, 'md', 'temp_equil', 'equilibrated.pdb')
    new_dcd_filename = 'umbrella1.dcd'
    new_pdb_filename = 'umbrella1.pdb'
    milestone.box_vectors = None
  else: # then some already exist
    number_list = []
    for existing_file in existing_umbrella_files:
      number_list.append(int(re.findall(r".+(\d+).dcd", existing_file)[0]))
    current_num = max(number_list)
    milestone.openmm.umbrella_pdb_filename = os.path.join(seekr_calc.project.rootdir, milestone.directory, 'md', 'umbrella', 'umbrella%d.pdb' % current_num)
    next_num = current_num + 1
    new_dcd_filename = 'umbrella%d.dcd' % next_num
    new_pdb_filename = 'umbrella%d.pdb' % next_num
  
  return new_dcd_filename, new_pdb_filename
    