'''
A sample script for running the forward and reverse stages

Created on May 15, 2018

@author: lvotapka
'''

import seekr
from seekr import amber
import sys, os
from simtk.unit import *
import mdtraj

print "Parse arguments"
which = None
if len(sys.argv) < 2: # then assume all
  which = 'all'
elif sys.argv[1] == 'all':
  which == 'all'
else:
  which = int(sys.argv[1])

print "Loading SEEKR calculation."

##################################################################
# VARIABLES WITHIN SECTION BELOW SHOULD BE MODIFIED TO YOUR SYSTEM
##################################################################


picklename = '/home/lvotapka/tryp_test/seekr_calc.pickle'
me = seekr.openSeekrCalc(picklename)

lig_selection = [3222, 3223, 3224, 3225, 3226, 3227, 3228, 3229, 3230]
rec_selection = [2466, 2478, 2489, 2535, 2745, 2769, 2787]

step_chunk_size = 1000
me.fwd_rev_stage.steps = step_chunk_size # in 2*fs
me.fwd_rev_stage.energy_freq = 1000
me.fwd_rev_stage.traj_freq = 1000
me.fwd_rev_stage.launches_per_config = 1
me.fwd_rev_stage.barostat = False # leave barostat off
transition_filename = 'transition_fwd.dat'

##################################################################
# DON'T MODIFY THE SECTION BELOW UNLESS YOU KNOW WHAT YOU'RE DOING
##################################################################

if which == 'all': # then run all milestones
  all_milestones = me.milestones
else:
  all_milestones = [me.milestones[which]]
  
for milestone in all_milestones:
  if milestone.md:
    if not milestone.openmm.prmtop_filename: 
      print "prmtop file not found for milestone %d. Skipping..." % milestone.index
      continue
      
    print "launching constant energy forward stage for milestone %d:" % which
    box_vectors = milestone.box_vectors
    milestone.atom_selection_1 = rec_selection
    milestone.atom_selection_2 = lig_selection
    fwd_rev_path = os.path.join(me.project.rootdir, milestone.directory, 'md', 'fwd_rev')
    
    me.fwd_rev_stage.success_coords_pickle = os.path.join(me.project.rootdir, milestone.directory, 'md', 'fwd_rev', 'success_coords.pickle') # TODO: remove line
    success_coords_pickle =  open(me.fwd_rev_stage.success_coords_pickle, 'rb') 
    success_coords_pickle_file = open(success_coords_pickle, 'rb')
    positions = pickle.load(success_coords_pickle_file)
    success_coords_pickle_file.close()
    
    me.fwd_rev_stage.success_vels_pickle = os.path.join(me.project.rootdir, me.milestones[which].directory, 'md', 'fwd_rev', 'success_vels.pickle') # TODO: remove line
    success_vels_pickle =  open(me.fwd_rev_stage.success_vels_pickle, 'rb') 
    success_vels_pickle_file = open(success_vels_pickle, 'rb')
    velocities = pickle.load(success_vels_pickle_file)
    success_vels_pickle_file.close()
    
    velocities = -1.0 * velocities # REVERSE velocities
    positions = iter(positions)
    
    traj_base = "forward"
    print "running forwards"
    success_positions, success_velocities, data_file_name, indices_list = seekr.launch_fwd_rev_stage(me, me.milestones[which], traj_base, False, positions, input_vels=velocities, box_vectors=box_vectors, transition_filename=transition_filename)
    
    
    # TODO: parse transition file information

