'''
Created on May 12, 2020

As OpenSEEKR performs its calculations, this module provides a number of tools
to analyze its transitions and trajectories.

Most notably, this module computes the rate constants and thermodynamics of
binding and unbinding for the ligand-receptor system.

@author: lvotapka
'''
import os
import glob

import numpy as np
import mdtraj

import seekr
import dig_deeper

def get_umbrella_avg_distance(me, milestone):
    umbrella_dir = os.path.join(me.project.rootdir, milestone.directory, 
                                'md', 'umbrella')
    prmtop = os.path.join(me.project.rootdir, milestone.directory, 'md', 
                          'building', 'holo.parm7')
    umbrella_dcd_list = glob.glob(os.path.join(umbrella_dir, 'umbrella*.dcd'))
    
    mytraj1 = mdtraj.load(umbrella_dcd_list, top=prmtop, 
                           atom_indices=milestone.atom_selection_1)
    mytraj2 = mdtraj.load(umbrella_dcd_list, top=prmtop, 
                           atom_indices=milestone.atom_selection_2)
    com1_array = mdtraj.compute_center_of_mass(mytraj1)
    com2_array = mdtraj.compute_center_of_mass(mytraj2)
    distance_total = 0.0
    for com1, com2 in zip(com1_array, com2_array):
        this_distance = np.linalg.norm(com2-com1)
        distance_total += this_distance
    
    avg_distance = 10.0 * distance_total / com1_array.shape[0]
    return avg_distance

def get_fwd_rev_avg_distance(me, milestone):
    fwd_rev_dir = os.path.join(me.project.rootdir, milestone.directory, 
                                'md', 'fwd_rev')
    prmtop = os.path.join(me.project.rootdir, milestone.directory, 'md', 
                          'building', 'holo.parm7')
    fwd_dcd_list = sorted(glob.glob(
        os.path.join(fwd_rev_dir, 'forward*.dcd')), 
        key=seekr.sort_forward_dcd_key)
    data_file_name = os.path.join(fwd_rev_dir, 'transition_fwd.dat')
    
    downward_indices = dig_deeper.read_data_file_transitions_down(
        data_file_name, destination='1')
    upward_indices = dig_deeper.read_data_file_transitions_down(
        data_file_name, destination='3')
    
    upward_distance_total = 0.0
    upward_counter = 0
    downward_distance_total = 0.0
    downward_counter = 0
    for i, dcd_file in enumerate(fwd_dcd_list):
        last_frame1 = seekr.load_last_mdtraj_frame(
            dcd_file, prmtop, atom_indices=milestone.atom_selection_1)
        last_frame2 = seekr.load_last_mdtraj_frame(
            dcd_file, prmtop, atom_indices=milestone.atom_selection_2)
        com1_array = mdtraj.compute_center_of_mass(last_frame1)
        com2_array = mdtraj.compute_center_of_mass(last_frame2)
        #print('com1_array:', com1_array, 'com2_array:', com2_array)
        this_distance = np.linalg.norm(com2_array[0]-com1_array[0])
        
        if i in upward_indices:
            upward_distance_total += this_distance
            upward_counter += 1
        elif i in downward_indices:
            downward_distance_total += this_distance
            downward_counter += 1
        if upward_counter >= 10 and downward_counter >= 10:
            # enough statistics
            break
    if upward_counter > 0:
        upward_avg_distance = 10.0 * upward_distance_total / upward_counter
    else:
        upward_avg_distance = None
        
    if downward_counter > 0:
        downward_avg_distance = 10.0 * downward_distance_total / \
            downward_counter
    else:
        downward_avg_distance = None
        
    return upward_avg_distance, downward_avg_distance
    
    