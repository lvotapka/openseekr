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