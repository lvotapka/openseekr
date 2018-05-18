'''
A sample script for running umbrella sampling

Created on May 15, 2018

@author: lvotapka
'''

import seekr
from seekr import amber
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import sys, os, time
from sys import stdout

verbose = True

def create_forces(seekrcalc, milestone, system):
  '''
  Add the umbrella force: which maintains the ligand on the surface of the 
  spherical milestone.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to run the simulation for
   - system: the OpenMM system object to add the force to
  Output:
  '''
  new_force = CustomCentroidBondForce(2, 'k*(distance(g1,g2)-r0)^2')
  k = new_force.addGlobalParameter('k', seekrcalc.umbrella_stage.force_constant)
  r0 = new_force.addGlobalParameter('r0', milestone.radius*angstrom)
  g1 = new_force.addGroup(milestone.atom_selection_1)
  g2 = new_force.addGroup(milestone.atom_selection_2)
  if verbose: print "k:", k, "r0:", r0, "g1:", g1, "g2:", g2
  new_force.addBond([g1, g2], [])
  if verbose: print "new_force.getNumGlobalParameters():", new_force.getNumGlobalParameters()
  if verbose: print "new_force.getNumPerBondParameters():", new_force.getNumPerBondParameters()
  system.addForce(new_force)
  return

def launch_umbrella_stage(seekrcalc, milestone, box_vectors=None, traj_name='umbrella1.dcd'):
  '''launch an umbrella sampling job.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to run the simulation for
  Output:
   - ending_box_vectors: An array of three vectors that represents the ending
   state of the periodic box. This may not be the same as was started, but could
   have changed through the course of a constant pressure simulation.
  '''
  prmtop_filename = milestone.openmm.prmtop_filename
  pdb_filename = milestone.openmm.umbrella_pdb_filename
  inpcrd_filename = milestone.openmm.inpcrd_filename
  if verbose: print "opening files:", prmtop_filename, inpcrd_filename
  prmtop = AmberPrmtopFile(prmtop_filename)
  inpcrd = AmberInpcrdFile(inpcrd_filename)
  pdb = PDBFile(pdb_filename)
  
  system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds) # This is fine because h-bonds are always constrained in water!
  integrator = LangevinIntegrator(me.master_temperature*kelvin, 1/picosecond, 0.002*picoseconds) #LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
  platform = Platform.getPlatformByName('CUDA')
  properties = seekrcalc.openmm.properties #{'CudaDeviceIndex':'0', 'CudaPrecision':'double'}
  
  # add restraints
  create_forces(seekrcalc, milestone, system) #system, prmtop.topology, inpcrd.positions, seekrcalc.min_equil.constrained)
  if seekrcalc.umbrella_stage.barostat:
    barostat = MonteCarloBarostat(seekrcalc.umbrella_stage.barostat_pressure, seekrcalc.master_temperature*kelvin, seekrcalc.umbrella_stage.barostat_freq)
    system.addForce(barostat)
  
  simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
  simulation.context.setPositions(pdb.positions)
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
  simulation.step(seekrcalc.umbrella_stage.steps)
  print "time:", time.time() - starttime, "s"
  end_state = simulation.context.getState(getPositions=True)
  ending_box_vectors = end_state.getPeriodicBoxVectors()
  milestone.openmm.simulation = simulation
  return ending_box_vectors

print "Parse arguments"
which = None
if len(sys.argv) < 2: # then assume all
  which = 'all'
elif sys.argv[1] == 'all':
  which == 'all'
else:
  which = int(sys.argv[1])

print "Loading SEEKR calculation."
rootdir = '/home/lvotapka/tryp_test'
num_NPT_steps = 10000
num_NVT_steps = 100000

picklename = os.path.join(rootdir, 'seekr_calc.pickle')
me = seekr.openSeekrCalc(picklename)
me.umbrella_stage.force_constant = 90.0*kilocalories_per_mole/angstroms**2
me.umbrella_stage.steps = num_NPT_steps # in 2*fs
me.umbrella_stage.energy_freq = 1000
me.umbrella_stage.traj_freq = 1000
me.umbrella_stage.barostat = True # turn on barostat, run in NPT
me.umbrella_stage.barostat_freq = 25
me.umbrella_stage.barostat_pressure = 1.0*bar

if which == 'all': # then run all milestones
  for milestone in me.milestones:
    if milestone.md:
      print "launching constant pressure umbrella sampling for milestone %d:" % which
      box_vectors = launch_umbrella_stage(me, milestone, traj_name='umbrella1.dcd')
      pdb_filename = os.path.join(me.project.rootdir, milestone.directory, 'md', 'umbrella', 'umbrella1.pdb')
      amber.save_restart(me, milestone, pdb_filename)
      milestone.openmm.umbrella_pdb_filename = pdb_filename
      me.umbrella_stage.barostat = False # turn off the barostat and run in NVT
      me.umbrella_stage.steps = num_NVT_steps # in 2*fs
      print "running the second stage of umbrella sampling."
      box_vectors2 = launch_umbrella_stage(me, milestone, box_vectors, 'umbrella2.dcd')
else:
  print "launching constant pressure umbrella sampling for milestone %d:" % which
  box_vectors = launch_umbrella_stage(me, me.milestones[which], traj_name='umbrella1.dcd')
  pdb_filename = os.path.join(me.project.rootdir, me.milestones[which].directory, 'md', 'umbrella', 'umbrella1.pdb')
  amber.save_restart(me, me.milestones[which], pdb_filename)
  me.milestones[which].openmm.umbrella_pdb_filename = pdb_filename
  me.umbrella_stage.barostat = False # turn off the barostat and run in NVT
  me.umbrella_stage.steps = num_NVT_steps # in 2*fs
  print "running the second stage of umbrella sampling."
  box_vectors2 = launch_umbrella_stage(me, me.milestones[which], box_vectors, 'umbrella2.dcd')
  
