'''
A sample script for running the reversal stage of a system

Created on May 15, 2018

@author: lvotapka
'''

import seekr
from seekr import amber
import sys, os, math
from simtk.unit import *
import mdtraj
from simtk.openmm.app import AmberInpcrdFile

print "Parse arguments"
which = None
if len(sys.argv) < 2: # then assume all
  which = 'all'
elif sys.argv[1] == 'all':
  which = 'all'
else:
  which = int(sys.argv[1])

if len(sys.argv) == 3:
  launches_per_config = int(sys.argv[2])
else:
  launches_per_config = 1
  
print "which:", which

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
umbrella_glob = 'umbrella*.dcd'
reversal_frames = (1002, 10002, 1)
pos_vel_chunk_size = 400
transition_filename = 'transition_rev.dat'
me.openmm.properties = {'CudaDeviceIndex':'1', 'CudaPrecision':'mixed'}

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
      prmtop_path = os.path.join(me.project.rootdir, milestone.directory, 'md', 'building', 'holo.parm7')
      inpcrd_path = os.path.join(me.project.rootdir, milestone.directory, 'md', 'building', 'holo.rst7')
      if os.path.exists(prmtop_path) and os.path.exists(inpcrd_path):
        milestone.openmm.prmtop_filename = prmtop_path
        milestone.openmm.inpcrd_filename = inpcrd_path
        inpcrd = AmberInpcrdFile(inpcrd_path)
        milestone.box_vectors = inpcrd.boxVectors
        print "box_vectors:", milestone.box_vectors
      else:
        print "prmtop or inpcrd file not found for milestone %d. Skipping..." % milestone.index
        continue
    print "launching constant energy reverse stage for milestone:", which
    box_vectors = milestone.box_vectors
    milestone.atom_selection_1 = rec_selection
    milestone.atom_selection_2 = lig_selection
    fwd_rev_path = os.path.join(me.project.rootdir, milestone.directory, 'md', 'fwd_rev')
    umbrella_traj = os.path.join(me.project.rootdir, milestone.directory, 'md', 'umbrella', umbrella_glob)
    parm_file_name = os.path.join(me.project.rootdir, milestone.directory, 'md', 'building', 'holo.parm7')
    trajout = os.path.join(me.project.rootdir, milestone.directory, 'md', 'umbrella', 'imaged.dcd')
    cpptraj_script_location = os.path.join(me.project.rootdir, milestone.directory, 'md', 'umbrella', 'image_umbrella.cpptraj')
    box_info = seekr.make_box_info(box_vectors)
    seekr.autoimage_traj(parm_file_name, umbrella_traj, trajout, box_info, cpptraj_script_location=cpptraj_script_location, writing_frames=reversal_frames)
    dcd = mdtraj.iterload(trajout, top=parm_file_name, chunk=1)
    traj_base = "reverse"
    print "running reversals"

    #num_frames = launches_per_config*(reversal_frames[1] - reversal_frames[0]) / reversal_frames[2]
    #print "num_frames:", num_frames
    #print "pos_vel_chunk_size:", pos_vel_chunk_size
    #print "math.ceil((1.0*num_frames) / pos_vel_chunk_size):", math.ceil((1.0*num_frames) / pos_vel_chunk_size)
    
    # TODO: PROBLEM! What do to about the existing transitions.dat file???
    
    #for i in range(int(math.ceil((1.0*num_frames) / pos_vel_chunk_size))):
    complete = False
    i = 0
    save_fwd_rev = False
    while not complete:
      print "Running chunk %d" % i
      success_positions, success_velocities, data_file_name, indices_list, complete = seekr.launch_fwd_rev_stage(me, milestone, traj_base, True, dcd, pos_vel_chunk_size, box_vectors=box_vectors, suffix='_%d' % i, save_fwd_rev=save_fwd_rev)
      save_fwd_rev = True
      if len(success_positions) == 0:
        print "Reversal stage failed for this chunk: No successful reversal trajectories completed."
      else:
        print "saving coordinates and velocities for the reversal stage. len(success_positions)", len(success_positions), "len(success_velocities):", len(success_velocities)
        seekr.pickle_coords_vels(me, milestone, success_positions, success_velocities, index=i)
      i += 1
    
me.save()
    
