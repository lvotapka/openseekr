'''
A sample script for running the forward and reverse stages

Created on May 15, 2018

@author: lvotapka
'''

import seekr
from seekr import amber
import sys, os, glob, re
from simtk.unit import *
import mdtraj
import pickle as pickle
from pprint import pprint
import numpy as np
from simtk.openmm.app import AmberInpcrdFile

print("Parse arguments")
which = None
if len(sys.argv) < 2: # then assume all
    which = 'all'
elif sys.argv[1] == 'all':
    which = 'all'
else:
    which = int(sys.argv[1])

print("Loading SEEKR calculation.")

##################################################################
# VARIABLES WITHIN SECTION BELOW SHOULD BE MODIFIED TO YOUR SYSTEM
##################################################################


picklename = '/home/lvotapka/tryp_test/seekr_calc.xml'
me = seekr.openSeekrCalc(picklename)

step_chunk_size = 1000
me.fwd_rev_stage.steps = step_chunk_size # in 2*fs
me.fwd_rev_stage.energy_freq = 1000
me.fwd_rev_stage.traj_freq = 1000
me.fwd_rev_stage.launches_per_config = 1
me.fwd_rev_stage.barostat = False # leave barostat off
transition_filename = 'transition_fwd.dat'
me.openmm.properties = {'CudaDeviceIndex':'0', 'CudaPrecision':'mixed'}

##################################################################
# DON'T MODIFY THE SECTION BELOW UNLESS YOU KNOW WHAT YOU'RE DOING
##################################################################

def add_dictionaries(dict1, dict2):
    '''
    adds the values numerically within each dictionary
    NOTE: dict1 is updated and returned BY REFERENCE
    '''
    new_dict = dict1
    for key in list(dict2.keys()):
        if key in list(dict1.keys()):
            dict1[key] += dict2[key]
        else:
            dict1[key] = dict2[key]

    return dict1

if which == 'all': # then run all milestones
    all_milestones = me.milestones
else:
    all_milestones = [me.milestones[which]]

for milestone in all_milestones:
    if milestone.md:
        #if not milestone.openmm.prmtop_filename:
        #  print "prmtop file not found for milestone %d. Skipping..." % milestone.index
        #  continue
        if not milestone.openmm.prmtop_filename:
            prmtop_path = os.path.join(me.project.rootdir, milestone.directory, 'md', 'building', 'holo.parm7')
            inpcrd_path = os.path.join(me.project.rootdir, milestone.directory, 'md', 'building', 'holo.rst7')
            if os.path.exists(prmtop_path) and os.path.exists(inpcrd_path):
                milestone.openmm.prmtop_filename = prmtop_path
                milestone.openmm.inpcrd_filename = inpcrd_path
                inpcrd = AmberInpcrdFile(inpcrd_path)
                milestone.fwd_rev_box_vectors = inpcrd.boxVectors
                print("box_vectors:", milestone.fwd_rev_box_vectors)
            else:
                print("prmtop or inpcrd file not found for milestone %d. Skipping..." % milestone.index)
                continue

        print("launching constant energy forward stage for milestone:", which)
        box_vectors = milestone.fwd_rev_box_vectors
        fwd_rev_path = os.path.join(me.project.rootdir, milestone.directory, 'md', 'fwd_rev')

        traj_base = "forward"
        transition_dict_total = {}
        incubation_time_list_total = []
        save_fwd_rev = False

        glob_list = sorted(glob.glob(os.path.join(me.project.rootdir,
                  milestone.directory, 'md', 'fwd_rev', 'success_coords*.pickle')),
                  key=seekr.sort_pickle_key)

        for coords_pickle in glob_list:
            i = int(re.search('\d+', os.path.basename(coords_pickle)).group(0))
            vels_pickle = os.path.join(me.project.rootdir, milestone.directory, 'md', 'fwd_rev', 'success_vels%d.pickle' % i)
            assert os.path.exists(vels_pickle), coords_pickle+' file exists, but corresponding velocity file: '+vels_pickle+' does not.'
            #me.fwd_rev_stage.success_coords_pickle = os.path.join(me.project.rootdir, milestone.directory, 'md', 'fwd_rev', 'success_coords.pickle') # TODO: remove line
            print("Opening pickles:", coords_pickle, vels_pickle)
            success_coords_pickle =  open(coords_pickle, 'rb')
            #success_coords_pickle_file = open(success_coords_pickle, 'rb')
            positions = pickle.load(success_coords_pickle)
            success_coords_pickle.close()

            #me.fwd_rev_stage.success_vels_pickle = os.path.join(me.project.rootdir, milestone.directory, 'md', 'fwd_rev', 'success_vels.pickle') # TODO: remove line
            success_vels_pickle =  open(vels_pickle, 'rb')
            #success_vels_pickle_file = open(success_vels_pickle, 'rb')
            velocities = pickle.load(success_vels_pickle)
            success_vels_pickle.close()

            reversed_vels = []
            print("Reversing Velocities")
            for vel in velocities:
                reversed_vels.append(-1.0 * vel)

            assert len(positions) == len(velocities), "The length of provided velocities and positions must be equal."
            positions = iter(positions)

            print("Running Forwards")
            complete = False

            while not complete:
                print("Running chunk %d" % i)

                success_positions, success_velocities, data_file_name, indices_list, complete = seekr.launch_fwd_rev_stage(me, milestone,
                              traj_base, False, positions, input_vels=reversed_vels, box_vectors=box_vectors,
                              transition_filename=transition_filename, suffix='_%d' % i, save_fwd_rev=save_fwd_rev)
                save_fwd_rev = True
                del success_positions
                del success_velocities
                del positions
                del velocities

                # TODO: parse transition file information
                transition_dict, avg_incubation_time, incubation_time_list = seekr.read_data_file_transitions(data_file_name, me, milestone)
                incubation_time_list_total += incubation_time_list
                transition_dict_total = add_dictionaries(transition_dict_total, transition_dict)
                print("TRANSITION DICT:")
                pprint(transition_dict)

        print("Transition Data:")
        print("transition dictionary:")
        pprint(transition_dict_total)
        print("average incubation time:")
        avg_incubation_time = np.mean(incubation_time_list_total)
        pprint(avg_incubation_time)

        seekr.pickle_transition_info(me, milestone, transition_dict_total, avg_incubation_time)
