'''
A sample script for running the reversal stage of a system

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


rootdir = '/home/lvotapka/tryp_test/seekr_calc.pickle'
me = seekr.openSeekrCalc(picklename)

lig_selection = [3222, 3223, 3224, 3225, 3226, 3227, 3228, 3229, 3230]
rec_selection = [2466, 2478, 2489, 2535, 2745, 2769, 2787]

step_chunk_size = 1000
me.fwd_rev_stage.steps = step_chunk_size # in 2*fs
me.fwd_rev_stage.energy_freq = 1000
me.fwd_rev_stage.traj_freq = 1000
me.fwd_rev_stage.launches_per_config = 10
me.fwd_rev_stage.barostat = False # leave barostat off
umbrella_glob = 'umbrella*.dcd'
reversal_frames = (1010, 10010, 1)

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
    print "launching constant energy reverse stage for milestone %d:" % which
    box_vectors = milestone.box_vectors
    milestone.atom_selection_1 = rec_selection
    milestone.atom_selection_2 = lig_selection
    fwd_rev_path = os.path.join(me.project.rootdir, milestone.directory, 'md', 'fwd_rev')
    umbrella_traj = os.path.join(me.project.rootdir, milestone.directory, 'md', 'umbrella', umbrella_glob)
    parm_file_name = os.path.join(me.project.rootdir, me.milestones[which].directory, 'md', 'building', 'holo.parm7')
    trajout = os.path.join(me.project.rootdir, me.milestones[which].directory, 'md', 'umbrella', 'imaged.dcd')
    cpptraj_script_location = os.path.join(me.project.rootdir, me.milestones[which].directory, 'md', 'umbrella', 'image_umbrella.cpptraj')
    box_info = seekr.make_box_info(box_vectors)
    seekr.autoimage_traj(parm_file_name, umbrella_traj, trajout, box_info, cpptraj_script_location=cpptraj_script_location, writing_frames=reversal_frames)
    
    dcd = mdtraj.load(trajout, top=parm_file_name)
    positions = dcd.xyz
    traj_base = "reverse"
    print "running reversals"
    starting_positions, starting_velocities, data_file_name, indices_list = seekr.launch_fwd_rev_stage(me, milestone, traj_base, True, positions, box_vectors=box_vectors)
    print "processing reversal data. len(starting_positions:", len(starting_positions), "len(starting_velocities):", len(starting_velocities)
    success_coords, success_vels = seekr.process_reversal_data(starting_positions, starting_velocities, data_file_name)
    print "saving coordinates and velocities for the forward stage. len(success_coords)", len(success_coords), "len(success_vels):", len(success_vels)
    seekr.pickle_coords_vels(me, milestone, starting_positions, starting_velocities, success_coords, success_vels)
    
    me.save()
    
  

