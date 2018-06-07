'''
A sample script for running the forward and reverse stages

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

from seekrplugin import SeekrForce

verbose = True

def create_seekr_force(seekrcalc, milestone, system):
  '''create the SEEKR 'force' even though it's more of a monitor than a force.'''
  force = SeekrForce()
  neighbor1 = seekrcalc.milestones[milestone.neighbors[0]]
  neighbor2 = seekrcalc.milestones[milestone.neighbors[1]]
  radius1 = neighbor1.radius
  radius2 = milestone.radius
  radius3 = neighbor2.radius
  print "mark0"
  print "milestone.atom_selection_1:", milestone.atom_selection_1
  print "milestone.atom_selection_2:", milestone.atom_selection_2
  force.addSphericalMilestone(len(milestone.atom_selection_1), len(milestone.atom_selection_2), radius1, radius2, radius3, milestone.atom_selection_1, milestone.atom_selection_2)
  print "force created"
  system.addForce(force)
  print "force added to system."
  

def launch_fwd_rev_stage(seekrcalc, milestone, box_vectors=None, traj_name='fwd_rev1.dcd'):
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
  
  # TODO: change this back
  properties = seekrcalc.openmm.properties #{'CudaDeviceIndex':'0', 'CudaPrecision':'mixed', 'UseCpuPme':'false'}
  
  # add restraints
  #create_forces(seekrcalc, milestone, system) #system, prmtop.topology, inpcrd.positions, seekrcalc.min_equil.constrained)
  create_seekr_force(seekrcalc, milestone, system)
  
  
  simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
  simulation.context.setPositions(pdb.positions)
  if box_vectors:
    simulation.context.setPeriodicBoxVectors(*box_vectors)
  elif inpcrd.boxVectors is not None:
    simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
  
  if verbose: print "Running energy minimization on milestone:", milestone.index
  simulation.minimizeEnergy() # TODO: remove this line once debugging is complete
  
  fwd_rev_traj = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'fwd_rev', traj_name)
  simulation.reporters.append(StateDataReporter(stdout, seekrcalc.fwd_rev_stage.energy_freq, step=True, potentialEnergy=True, temperature=True, volume=True))
  simulation.reporters.append(DCDReporter(fwd_rev_traj, seekrcalc.fwd_rev_stage.traj_freq))
  starttime = time.time()
  simulation.step(seekrcalc.fwd_rev_stage.steps)
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
num_steps = 1000

picklename = os.path.join(rootdir, 'seekr_calc.pickle')
me = seekr.openSeekrCalc(picklename)
me.fwd_rev_stage.steps = num_steps # in 2*fs
me.fwd_rev_stage.energy_freq = 1000
me.fwd_rev_stage.traj_freq = 1000
me.fwd_rev_stage.barostat = False # turn on barostat, run in NPT

if which == 'all': # then run all milestones
  raise Exception, "Enter a number as argument for now."
  for milestone in me.milestones:
    if milestone.md:
      print "launching constant energy forward/reverse stage for milestone %d:" % which
      box_vectors = launch_fwd_rev_stage(me, milestone, traj_name='fwd_rev1.dcd')
      pdb_filename = os.path.join(me.project.rootdir, milestone.directory, 'md', 'fwd_rev', 'fwd_rev1.pdb')
      amber.save_restart(me, milestone, pdb_filename)
      milestone.openmm.fwd_rev_pdb_filename = pdb_filename
      me.fwd_rev_stage.steps = num_steps # in 2*fs
      #print "running the second stage of umbrella sampling."
      #box_vectors2 = launch_fwd_rev_stage(me, milestone, box_vectors, 'fwd_rev2.dcd')
else:
  print "launching constant energy forward/reverse stage for milestone %d:" % which
  box_vectors = launch_fwd_rev_stage(me, me.milestones[which], traj_name='fwd_rev1.dcd')
  '''
  pdb_filename = os.path.join(me.project.rootdir, me.milestones[which].directory, 'md', 'fwd_rev', 'fwd_rev1.pdb')
  amber.save_restart(me, me.milestones[which], pdb_filename)
  me.milestones[which].openmm.fwd_rev_pdb_filename = pdb_filename
  me.fwd_rev_stage.steps = num_steps # in 2*fs
  
  #box_vectors2 = launch_umbrella_stage(me, me.milestones[which], box_vectors, 'umbrella2.dcd')
  '''
  


