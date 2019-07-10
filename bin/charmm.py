'''
Merge ligand and receptor files to write parameter/topology/coordinate files
for Charmm simulations.

Created on Jun 11, 2019

@author: lvotapka
'''

import os, string
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from copy import deepcopy
import parmed
import psfmerge
import amber

verbose = True

class CharmmSettings():
  def __init__(self):
    self.receptor_psf = ''
    self.ligand_psf = ''
    self.param_list = []
    self.insert_index = -1

def charmm_building(seekrcalc, milestone, charmm_settings):
  '''the pre-minimization procedure for charmm ff simulations
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to prepare the AMBER system for
   ###- amber_settings: the AmberSettings() object that defines settings for the
       ###tleap calculation.
  Output:
   - None
  '''
  i = milestone.index
  building = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'building')
  pdb_filename = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'holo_wet.pdb')
  psf_filename = os.path.join(building, 'holo_wet.psf')
  
  if not os.path.exists(pdb_filename):
    print "Skipping milestone %d: no holo structure detected at %s." % (i, pdb_filename)
    return
  
  insert_index = charmm_settings.insert_index
  if os.path.exists(psf_filename):
    print "Skipping psfmerge: file already exists:", psf_filename
  else:
    newpsf = psfmerge.merge_psf_files(charmm_settings.receptor_psf, 
                                    charmm_settings.ligand_psf, 
                                    insert_index = insert_index) # merge the psf files
    newpsf.write(psf_filename) # write the new psf files
  
  milestone.openmm.psf_filename = psf_filename
  milestone.openmm.wet_holo_pdb_filename = pdb_filename
  milestone.openmm.charmm_params_filename_list = charmm_settings.param_list
  return

def create_simulation(seekrcalc, milestone, box_vectors, box_lengths_angles):
  '''create the OpenMM simulation object.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to prepare the CHARMM system for
  Output:
   - None
  '''
  if not milestone.openmm.psf_filename or not milestone.openmm.wet_holo_pdb_filename:
    print "Charmm parameter/topology/coordinate files not found for milestone: %d. Skipping." % milestone.index
    return
  psf_filename = milestone.openmm.psf_filename
  pdb_filename = milestone.openmm.wet_holo_pdb_filename
  #rtf_filename = milestone.openmm.rtf_filename
  #par_filename = milestone.openmm.par_filename
  all_param_list = milestone.openmm.charmm_params_filename_list
  
  mypdb = PDBFile(pdb_filename)
  psf = CharmmPsfFile(psf_filename)
  params = CharmmParameterSet(*all_param_list)
  #psf.setBox(6.183293595232648, 6.183293595232648, 9.477549753348887)
  psf.setBox(*box_lengths_angles)

  system = psf.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=1*nanometer)
  
  integrator = deepcopy(seekrcalc.min_equil.temp_equil_integrator) #LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
  platform = seekrcalc.openmm.platform #Platform.getPlatformByName('CUDA')
  properties = seekrcalc.openmm.properties #{'CudaDeviceIndex':'0', 'CudaPrecision':'double'}
  
  # add restraints
  amber.create_restraints(system, psf.topology, seekrcalc.min_equil.constrained) # TODO: clean this up
  #barostat = MonteCarloBarostat(1.0*bar, seekrcalc.master_temperature*kelvin, 25)
  barostat = MonteCarloMembraneBarostat(1.0*bar, 0.0*bar*nanometers, seekrcalc.master_temperature*kelvin, MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFree , 25)
  system.addForce(barostat)
  
  simulation = Simulation(psf.topology, system, integrator, platform, properties)
  simulation.context.setPositions(mypdb.positions)
  print "box_vectors:", box_vectors
  print "box_vectors[0]:", box_vectors[0]
  print "box_vectors[1]:", box_vectors[1]
  print "box_vectors[2]:", box_vectors[2]
  #simulation.context.setPeriodicBoxVectors(*box_vectors)
  simulation.context.setPeriodicBoxVectors(box_vectors[0], box_vectors[1], box_vectors[2])
  
  milestone.openmm.system = system
  milestone.openmm.simulation = simulation
  return
  
def save_restart(seekrcalc, milestone, pdb_save_filename=None):
  '''Save CHARMM files for easy restart.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to prepare the AMBER system for
  Output:
   - None
  '''
  if not milestone.openmm.psf_filename or not milestone.openmm.wet_holo_pdb_filename:
    print "Charmm parameter/topology/coordinate files not found for milestone: %d. Skipping." % milestone.index
    return
  state = milestone.openmm.simulation.context.getState(getPositions = True, enforcePeriodicBox = True)
  positions = state.getPositions()
  charmm_structure = parmed.load_file(milestone.openmm.wet_holo_pdb_filename)
  charmm_structure.positions = positions
  if not pdb_save_filename:
    pdb_save_filename = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'temp_equil', 'equilibrated.pdb')
  charmm_structure.save(pdb_save_filename, overwrite=True)
  return pdb_save_filename
  
