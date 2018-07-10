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
      
    print "launching constant energy forward stage for milestone:", which
    box_vectors = milestone.box_vectors
    milestone.atom_selection_1 = rec_selection
    milestone.atom_selection_2 = lig_selection
    fwd_rev_path = os.path.join(me.project.rootdir, milestone.directory, 'md', 'fwd_rev')
    
    me.fwd_rev_stage.success_coords_pickle = os.path.join(me.project.rootdir, milestone.directory, 'md', 'fwd_rev', 'success_coords.pickle') # TODO: remove line
    success_coords_pickle =  open(me.fwd_rev_stage.success_coords_pickle, 'rb') 
    #success_coords_pickle_file = open(success_coords_pickle, 'rb')
    positions = pickle.load(success_coords_pickle)
    success_coords_pickle.close()
    
    me.fwd_rev_stage.success_vels_pickle = os.path.join(me.project.rootdir, me.milestones[which].directory, 'md', 'fwd_rev', 'success_vels.pickle') # TODO: remove line
    success_vels_pickle =  open(me.fwd_rev_stage.success_vels_pickle, 'rb') 
    #success_vels_pickle_file = open(success_vels_pickle, 'rb')
    velocities = pickle.load(success_vels_pickle)
    success_vels_pickle.close()
    
    reversed_vels = []
    print "Reversing Velocities"
    for vel in velocities:
      reversed_vels.append(-1.0 * vel)
    
    assert len(positions) == len(velocities), "The length of provided velocities and positions must be equal."
    positions = iter(positions)
    
    traj_base = "forward"
    print "Running Forwards"
    success_positions, success_velocities, data_file_name, indices_list = seekr.launch_fwd_rev_stage(me, me.milestones[which], traj_base, False, positions, input_vels=reversed_vels, box_vectors=box_vectors, transition_filename=transition_filename)
    
    # TODO: parse transition file information
    transition_dict, avg_incubation_time = seekr.read_data_file_transitions(data_file_name, me, milestone)
    
    print "Transition Data:"
    print "transition dictionary:"
    pprint(transition_dict)
    print "average incubation time:"
    pprint(avg_incubation_time)
    
    seekr.pickle_transition_info(me, milestone, transition_dict, avg_incubation_time)
    
