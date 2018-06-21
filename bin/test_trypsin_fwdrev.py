'''
A sample script for running the forward and reverse stages

Created on May 15, 2018

@author: lvotapka
'''

import seekr
from seekr import amber
import sys, os
from simtk.unit import *

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



lig_selection = [3222, 3223, 3224, 3225, 3226, 3227, 3228, 3229, 3230]
rec_selection = [2466, 2478, 2489, 2535, 2745, 2769, 2787]

picklename = os.path.join(rootdir, 'seekr_calc.pickle')
me = seekr.openSeekrCalc(picklename)
step_chunk_size = 1000
me.fwd_rev_stage.steps = step_chunk_size # in 2*fs
me.fwd_rev_stage.energy_freq = 1000
me.fwd_rev_stage.traj_freq = 1000
me.fwd_rev_stage.launches_per_config = 10
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
  box_vectors = me.milestones[which].box_vectors
  me.milestones[which].atom_selection_1 = rec_selection
  me.milestones[which].atom_selection_2 = lig_selection
  umbrella_path = os.path.join(me.project.rootdir, me.milestones[which].directory, 'md', 'umbrella')
  trajin = os.path.join(umbrella_path, 'umbrella2.dcd')
  trajout = os.path.join(umbrella_path, 'umbrella2_imaged.dcd')
  autoimage_traj(me.milestones[which].openmm.prmtop_filename, trajin, trajout, make_box_info(box_vectors))
  dcd = mdtraj.load(trajout, top=os.path.join(me.project.rootdir, me.milestones[which].directory, 'md', 'building', 'holo.parm7'))
  positions = dcd.xyz[:10]
  box_vectors = launch_fwd_rev_stage(me, me.milestones[which], positions, box_vectors)
  '''
  pdb_filename = os.path.join(me.project.rootdir, me.milestones[which].directory, 'md', 'fwd_rev', 'fwd_rev1.pdb')
  amber.save_restart(me, me.milestones[which], pdb_filename)
  me.milestones[which].openmm.fwd_rev_pdb_filename = pdb_filename
  me.fwd_rev_stage.steps = num_steps # in 2*fs
  
  #box_vectors2 = launch_umbrella_stage(me, me.milestones[which], box_vectors, 'umbrella2.dcd')
  '''
  


