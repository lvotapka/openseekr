'''
This script extracts a successful forward stage trajectory that ends on a lower
milestone. It extracts the last frame of that trajectory and writes the
structure to the next milestone down's directory.

Created on June 28, 2018

@author: lvotapka
'''
import sys
import os
import glob
import pprint
from shutil import copyfile
from math import sqrt, pi, exp, log, sin, cos

import mdtraj
import numpy as np

import seekr

def get_com(coords):
    '''
    TODO: add docstring
    '''
    com = np.zeros([1, 1, 3])
    n = len(coords[0])
    for i in range(n):
        coord = coords[0, i, :]
        com += coord
    com = com / n
    return com

def get_com_offset(coords0, coords1):
    '''
    TODO: add docstring
    '''
    assert len(coords0[0]) == len(coords1[0]), "there needs to be an equal number of atoms selected in coords0 and coords1."
    com0 = np.zeros([1,1,3])
    com1 = np.zeros([1,1,3])
    n = len(coords0[0])
    for i in range(n):
        coord0 = coords0[0,i,:]
        com0 += coord0
        coord1 = coords1[0,i,:]
        com1 += coord1
    com0 = com0 / n
    com1 = com1 / n
    offset = com1 - com0
    return offset

def compute_rmsd(coords0, coords1, offset):
    '''
    #TODO: add docstring
    '''
    assert len(coords0[0]) == len(coords1[0]), "there needs to be an equal number of atoms selected in coords0 and coords1."
    sd = 0.0 # square deviation
    n = len(coords0[0])
    for i in range(n):
        coord0 = coords0[0,i,:] + offset
        coord1 = coords1[0,i,:]
        diff_vec = coord1 - coord0
        diff_vec_sq = float(diff_vec[0,0,0]**2 + diff_vec[0,0,1]**2 + diff_vec[0,0,2]**2)
        sd += diff_vec_sq

    msd = sd / n
    rmsd = sqrt(msd)
    return rmsd

def find_closest_ligand_orientation(prmtop, dcd_list, reference_top, reference_crd, ligname, center_atom_indices):
    '''Reads in the dcd_list files one by one, extracting their final frames, then
    compare them to a reference structure. The protein structures are aligned,
    then the center of masses of the dcd last frame and the reference are
    superimposed. Then the RMSD is measured. The frame with the smallest RMSD is
    extracted.
    Input:
     - prmtop: the prmtop file name
     - dcd_list: a list of string representing the locations of the dcd files
     - reference: the reference PDB file for this system
     - ligname: the RESNAME field in the PDB for the ligand
    Output:
     - best_last_frame: the last frame whose RMSD is closest to the reference
       structure
    '''
    best_last_frame = None
    best_rmsd = 9e9
    print("running...")
    ref_struct_all = mdtraj.load(reference_crd, top=reference_top)
    ref_top_all = ref_struct_all.topology
    query = "name CA or (resname '%s')" % ligname
    ref_alpha_carbons_and_lig = ref_top_all.select(query) # select protein alpha carbons

    last_frame_all = seekr.load_last_mdtraj_frame(dcd_list[0], prmtop)
    last_frame_top_all = last_frame_all.topology
    query2 = "(name CA or resname '%s')" % ligname
    last_frame_alpha_carbons_and_lig = last_frame_top_all.select(query2)

    ref_struct = mdtraj.load(reference_crd, top=reference_top, atom_indices=ref_alpha_carbons_and_lig)
    ref_top = ref_struct.topology
    ref_alpha_carbons = ref_top.select("name CA") # select protein alpha carbons
    ref_lig = ref_top.select("resname '%s'" % ligname)
    ref_coords = ref_struct.xyz[:,ref_lig]
    assert len(ref_lig) > 0, "no reference atoms selected for ligand with ligname %s" % ligname

    counter = 0
    for dcd_file in dcd_list:
        last_frame = seekr.load_last_mdtraj_frame(dcd_file, prmtop, atom_indices=last_frame_alpha_carbons_and_lig)
        last_frame_site = seekr.load_last_mdtraj_frame(dcd_file, prmtop, atom_indices=[x-1 for x in center_atom_indices])
        last_frame_top = last_frame.topology
        last_frame_site_top = last_frame_site.topology

        last_frame_alpha_carbons = last_frame_top.select('name CA')
        last_frame_lig = last_frame_top.select("resname '%s'" % ligname)
        last_frame_coords = last_frame.xyz[:,last_frame_lig]
        #print "last_frame_coords:", last_frame_coords
        assert len(ref_lig) > 0, "no last_frame atoms selected for ligand with ligname %s" % ligname
        assert len(ref_alpha_carbons) == len(last_frame_alpha_carbons), 'The reference and the dcd must have the same number of alpha carbons'
        assert len(ref_lig) == len(last_frame_lig), 'The reference and the dcd must have the same number of ligand atoms'

        last_frame_site_coords = last_frame_site.xyz[:]
        site_com = get_com(last_frame_site_coords)
        lig_com = get_com(last_frame_coords)
        dist_from_origin = np.sqrt((lig_com[0, 0, 0]-site_com[0, 0, 0])**2 + (lig_com[0, 0, 1]-site_com[0, 0, 1])**2 + (lig_com[0, 0, 2]-site_com[0, 0, 2])**2)

        last_frame.superpose(reference=ref_struct, frame=0, atom_indices=last_frame_alpha_carbons, ref_atom_indices=ref_alpha_carbons)
        #rmsd=mdtraj.rmsd(target=last_frame, reference=ref_struct, frame=0, atom_indices=last_frame_lig, ref_atom_indices=ref_lig)
        last_frame_coords = last_frame.xyz[:,last_frame_lig]

        offset = get_com_offset(last_frame_coords, ref_coords)
        rmsd = compute_rmsd(last_frame_coords, ref_coords, offset)

        #print "len(ref_coords):", len(ref_coords[0])
        #print "ref_coords:", ref_coords
        #print "len(last_frame_coords):", len(last_frame_coords[0])
        #print "last_frame_coords:", last_frame_coords

        print("rmsd of frame %i" % counter, rmsd, "distance from milestone origin:", dist_from_origin, "nm")
        if rmsd < best_rmsd:
            best_rmsd = rmsd
            best_counter = counter
            best_last_frame_alphas = last_frame

        counter += 1

        #ref_struct.save_pdb('/tmp/refframe.pdb') # DEBUG lines
        #last_frame.save_pdb('/tmp/lastframe.pdb') # DEBUG lines

    best_last_frame = seekr.load_last_mdtraj_frame(dcd_list[best_counter], prmtop)

    print("best_rmsd:", best_rmsd)

    return best_last_frame


def read_data_file_transitions_down(data_file_name, destination='1', last_frame=True):
    '''Read transition data file, return the first instance of a transition to a
    lower milestone. Return the file index of the forward trajectory to this lower
    milestone.
    Input:
     - data_file_name: string of the filename of a milestone's forward transition
     file
     - destination: the string to search for. In the downward milestone case: '1'.
        Upward would be '3'.
    Output:
     - downward_index: integer of the forward trajectory that ends on the lower
     milestone
    '''
    downward_indices = []
    data_file = open(data_file_name, 'r')
    for i, line in enumerate(data_file.readlines()):
        line = line.split()
        #print "line:", line, "destination:", destination
        if line[0] == destination:
            downward_indices.append(i)
    data_file.close()
    assert len(downward_indices) != 0, "FAILURE: no downward forward trajectories detected. Please run additional umbrella sampling and rev/fwd trajectories."
    return downward_indices

if __name__ == "__main__":
    downward = True
    print("Parse arguments")
    if len(sys.argv) < 4:
        print("Usage:\npython dig_deeper.py MILESTONE SEEKRXML METHOD *ARGUMENTS")
        print("Available arguments for 'METHOD': first, last, similar, index")
        print("Usage for 'similar' method:")
        print("python dig_deeper.py MILESTONE SEEKRXML similar REF_PARM7 REF_RST7  LIG_RESNAME")
        print("Usage for 'index' method:")
        print("python dig_deeper.py MILESTONE SEEKRXML index FILENAME")
        print("be sure to provide reference PDB and ligand resname if using 'similar' method argument.")
        print("if the last argument is 'up' then an upward-going trajectory is chosen.")
        exit()
    
    which = int(sys.argv[1])
    picklename = sys.argv[2]
    method = sys.argv[3]
    ref_pdb = None
    lig_resname = None
    
    if method == 'similar':
        #ref_pdb = sys.argv[4]
        ref_parm7 = sys.argv[4]
        ref_rst7 = sys.argv[5]
        lig_resname = sys.argv[6]
    elif method == 'index':
        index = sys.argv[4]
    
    if sys.argv[-1] == 'up':
        downward = False
    
    print("Loading SEEKR calculation.")
    me = seekr.openSeekrCalc(picklename)
    
    milestone = me.milestones[which]
    if downward:
        lower_milestone = me.milestones[which-1] # TODO: hacky
    else:
        lower_milestone = me.milestones[which+1]
    
    # define all directories and files
    fwd_rev_dir = os.path.join(me.project.rootdir, milestone.directory, 'md', 'fwd_rev')
    lower_temp_equil_dir = os.path.join(me.project.rootdir, lower_milestone.directory, 'md', 'temp_equil')
    lower_temp_equil_filename = os.path.join(lower_temp_equil_dir, 'equilibrated.pdb')
    lower_milestone_holo = os.path.join(me.project.rootdir, lower_milestone.directory, 'md', 'holo_wet.pdb')
    lower_milestone_building = os.path.join(me.project.rootdir, lower_milestone.directory, 'md', 'building')
    data_file_name = os.path.join(fwd_rev_dir, 'transition_fwd.dat') # find the index of a successful downward trajectory
    prmtop = os.path.join(me.project.rootdir, milestone.directory, 'md', 'building', 'holo.parm7')
    new_prmtop = os.path.join(lower_milestone_building, 'holo.parm7')
    new_inpcrd = os.path.join(lower_milestone_building, 'holo.rst7')
    
    # NEXT TIME: find the downward forward dcd filenames for 'similar' method
    
    # figure out which forward to pull out from
    print("Attempting to extract the last frame of a successful downward trajectory.")
    if downward:
        downward_indices = read_data_file_transitions_down(data_file_name)
    else:
        downward_indices = read_data_file_transitions_down(data_file_name, destination='3') # TODO: hacky
    
    dcd_list = sorted(glob.glob(os.path.join(fwd_rev_dir, 'forward*.dcd')), key=seekr.sort_forward_dcd_key)
    
    print("Writing new structures and files needed to run umbrella simulation on the lower milestone (milestone %d)" % lower_milestone.index)
    #last_fwd_frame = mdtraj.load(downward_fwd_dcd, top=prmtop)[-1]
    if method=='first':
        #downward_fwd_dcd = os.path.join(fwd_rev_dir, 'forward%i_0.dcd' % downward_indices[0])
        downward_dcd = dcd_list[downward_indices[0]]
        print("Extracting frame from file:", downward_dcd)
        last_fwd_frame = seekr.load_last_mdtraj_frame(downward_dcd, prmtop)
    elif method=='last':
        #downward_fwd_dcd = os.path.join(fwd_rev_dir, 'forward%i_0.dcd' % downward_indices[-1])
        downward_dcd = dcd_list[downward_indices[-1]]
        print("Extracting frame from file:", downward_dcd)
        last_fwd_frame = seekr.load_last_mdtraj_frame(downward_dcd, prmtop)
    elif method=='similar':
        dcd_downward_list = []
        for dcd_index in downward_indices:
            #dcd_list.append(os.path.join(fwd_rev_dir, 'forward%i_0.dcd' % dcd_index))
            dcd_downward_list.append(dcd_list[dcd_index])
        last_fwd_frame = find_closest_ligand_orientation(prmtop, dcd_downward_list, ref_parm7, ref_rst7, lig_resname, lower_milestone.center_atom_indices)
    elif method=='index':
        downward_dcd = os.path.join(fwd_rev_dir, index)
        print("Extracting frame from file:", downward_dcd)
        last_fwd_frame = seekr.load_last_mdtraj_frame(downward_dcd, prmtop)
    else:
        raise Exception("Method not allowed: %s" % method)
    
    last_fwd_frame.save_pdb(lower_temp_equil_filename)
    last_fwd_frame.save_pdb(lower_milestone_holo)
    last_fwd_frame.save_amberrst7(new_inpcrd)
    
    # copy the prmtop to the lower building directory
    copyfile(prmtop, new_prmtop)
    lower_milestone.openmm.prmtop_filename = new_prmtop
    lower_milestone.openmm.inpcrd_filename = new_inpcrd
    
    lower_milestone.temp_equil_box_vectors = milestone.fwd_rev_box_vectors
    
    me.save()
