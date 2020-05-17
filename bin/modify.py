'''
Created on Mar 15, 2019

@author: lvotapka

Module for modifying, adding, or deleting milestones
'''

import sys
import os
import math
import shutil
import subprocess
import glob

import numpy as np

import seekr
#from seekr import amber
from seekr import deserialize_transition_info
import analyze

quick_analysis = True

def add_milestone(me, radius):
    '''
    Add a new Concentric spherical milestone object to this seekrCalculation.
    '''
    siteid = 0
    spheres_in_site = 0
    insert_index = None # index to insert before
    for i, milestone in enumerate(me.milestones):
        if milestone.radius > radius and insert_index == None:
            insert_index = i
            absolute = milestone.absolute
            atom_indices = milestone.center_atom_indices
            atom_selection_1 = milestone.atom_selection_1
            atom_selection_2 = milestone.atom_selection_2
            center_vec = milestone.center_vec

        if milestone.siteid == siteid:
            spheres_in_site += 1

    assert insert_index != None, 'The provided radius exists beyond the '\
        'current BD milestone. This is not currently supported. Move the '\
        'existing BD milestone, and add a new milestone below it.'

    # create the new milestone
    new_milestone = seekr.Concentric_Spherical_Milestone(i, siteid, absolute, 
                                                         md=me.project.md, 
                                                         bd=me.project.bd)
    anchor = np.array([0.0, 0.0, 0.0]) # TODO: may modification in the future
    new_milestone.index = insert_index
    new_milestone.center_atom_indices = atom_indices
    new_milestone.center_vec = center_vec
    new_milestone.radius = radius
    new_milestone.anchor = anchor
    new_milestone.atom_selection_1 = atom_selection_1
    new_milestone.atom_selection_2 = atom_selection_2
    lower_neighbor_index = None
    upper_neighbor_index = None
    if insert_index > 0:
        lower_neighbor_index = insert_index - 1
        lower_neighbor = me.milestones[lower_neighbor_index]
        new_milestone.neighbors.append(lower_neighbor)
        # Lower neighbor needs its neighbor modified
        me.milestones[insert_index-1].neighbors[1] = new_milestone
        me.milestones[insert_index].neighbors[0] = new_milestone
    elif insert_index == 0:
        me.milestones[insert_index].neighbors.insert(0, new_milestone)

    #upper_neighbor_index = insert_index + 1
    upper_neighbor = me.milestones[insert_index]
    new_milestone.neighbors.append(upper_neighbor)
    new_milestone.bd = False
    new_milestone.fullname = "%d_%d_%d_%.1f_%.1f_%.1f" % (insert_index+\
            (spheres_in_site*int(siteid)), insert_index, siteid, anchor[0], 
            anchor[1], anchor[2])
    new_milestone.directory = "anchor_%s" % new_milestone.fullname
    
    for i in range(insert_index, len(me.milestones)): # modify upper milestones
        me.milestones[i].index = i+1
        old_dirname = me.milestones[i].directory
        id1, id2, mysiteid, anchor0, anchor1, anchor2 = \
                me.milestones[i].fullname.split('_')
        new_fullname = "%d_%d_%s_%s_%s_%s" % (int(id1)+1, int(id2)+1, 
                                              mysiteid, anchor0, anchor1, 
                                              anchor2)
        new_dirname = 'anchor_%s' % new_fullname
        old_path = os.path.join(me.project.rootdir, old_dirname)
        new_path = os.path.join(me.project.rootdir, new_dirname)
        if me.milestones[i].openmm.prmtop_filename:
            me.milestones[i].openmm.prmtop_filename = os.path.join(new_path, 
                    'md', 'building','holo.parm7')
        if me.milestones[i].openmm.inpcrd_filename:
            me.milestones[i].openmm.inpcrd_filename = os.path.join(new_path, 
                    'md', 'building','holo.rst7')
            
    # insert new milestone into the milestone list in order
    me.milestones.insert(insert_index, new_milestone)
    return new_milestone, insert_index
    
def add_milestone_directory(me, new_milestone, insert_index):
    # modify neighbor milestones
    for i in range(insert_index+1, len(me.milestones)): 
        # modify upper milestones
        old_dirname = me.milestones[i].directory
        id1, id2, mysiteid, anchor0, anchor1, anchor2 = \
                me.milestones[i].fullname.split('_')
        new_fullname = "%d_%d_%s_%s_%s_%s" % (int(id1)+1, int(id2)+1, 
                mysiteid, anchor0, anchor1, anchor2)
        new_dirname = 'anchor_%s' % new_fullname
        old_path = os.path.join(me.project.rootdir, old_dirname)
        old_path_tmp = old_path+"temp"
        new_path = os.path.join(me.project.rootdir, new_dirname)
        new_path_tmp = new_path+"temp"
        
        me.milestones[i].fullname = new_fullname
        me.milestones[i].directory = new_dirname
        
        if os.path.exists(new_path):
            #tmp_name = os.path.join(me.project.rootdir, 'temp_dir')
            if os.path.exists(old_path_tmp):
                print("temporarily renaming", old_path_tmp, "to", new_path_tmp)
                os.rename(old_path_tmp, new_path_tmp)
            else:
                print("temporarily renaming", old_path, "to", new_path_tmp)
                os.rename(old_path, new_path_tmp)
        else:
            if os.path.exists(old_path_tmp):
                print("renaming files:", old_path_tmp, "to", new_path)
                os.rename(old_path_tmp, new_path)
            else:
                print("renaming files:", old_path, "to", new_path)
                os.rename(old_path, new_path)
    
    anchor_dir = os.path.join(me.project.rootdir, new_milestone.directory)
    anchor_filetree = seekr.Filetree({new_milestone.directory:{}})
    print("Making directory:", anchor_dir)
    anchor_filetree.make_tree(me.project.rootdir) # create anchor directory
    md_dir=os.path.join(anchor_dir,'md') # directory for MD
    md_filetree=seekr.Filetree({'md':seekr.mdtree})
    md_filetree.make_tree(anchor_dir) # make the MD filetree
    print("Making directory:", md_dir)
    return
    
def delete_milestone(me, index):
    delete_index = index
    assert delete_index >=0, "Delete index must be positive integer"
    assert delete_index < len(me.milestones)-1, "The provided delete index is "\
            "the BD milestone. This is not currently supported. Delete "\
            "interior milestone and move BD milestone to desired place."
            
    
    for i in range(delete_index+1, len(me.milestones)): 
        # modify upper milestones
        me.milestones[i].index = i-1
    
    del_milestone = me.milestones.pop(delete_index)
    if delete_index == 0:
        upper_neighbor = del_milestone.neighbors[0]
        upper_neighbor.neighbors.pop(0)
    else:
        upper_neighbor = del_milestone.neighbors[1]
        lower_neighbor = del_milestone.neighbors[0]
        upper_neighbor.neighbors[0] = lower_neighbor
        if lower_neighbor.index == 0:
            lower_neighbor.neighbors[0] = upper_neighbor
        else:
            lower_neighbor.neighbors[1] = upper_neighbor
    print("deleting milestone index:", del_milestone.index, "radius:", 
          del_milestone.radius)
    del_dir = os.path.join(me.project.rootdir, del_milestone.directory)
    return del_dir, delete_index

def delete_milestone_directory(del_dir, delete_index):
    shutil.rmtree(del_dir)
    for i in range(delete_index, len(me.milestones)): # modify upper milestones
        old_fullname = me.milestones[i].fullname
        old_dirname = me.milestones[i].directory
        id1, id2, mysiteid, anchor0, anchor1, anchor2 = \
                me.milestones[i].fullname.split('_')
        new_fullname = "%d_%d_%s_%s_%s_%s" % (int(id1)-1, int(id2)-1, 
                                              mysiteid, anchor0, anchor1, 
                                              anchor2)
        new_dirname = 'anchor_%s' % new_fullname
        #print "old_dirname:", old_dirname, "new_dirname:", new_dirname
        me.milestones[i].fullname = new_fullname
        me.milestones[i].directory = new_dirname
        old_path = os.path.join(me.project.rootdir, old_dirname)
        new_path = os.path.join(me.project.rootdir, new_dirname)
        print("renaming", old_path, "to", new_path)
        if me.milestones[i].openmm.prmtop_filename:
            me.milestones[i].openmm.prmtop_filename = os.path.join(new_path, 
                    'md', 'building','holo.parm7')
        if me.milestones[i].openmm.inpcrd_filename:
            me.milestones[i].openmm.inpcrd_filename = os.path.join(new_path, 
                    'md', 'building','holo.rst7')
        os.rename(old_path, new_path)
    
    return
    
def modify_milestone(me, index, radius):
    modify_index = index
    modify_radius = radius
    if modify_index == 0: # modifying the lowest milestone
        upper_neighbor = me.milestones[modify_index].neighbors[0]
        assert radius < upper_neighbor.radius, "Cannot choose a radius "\
            "larger than the upper adjacent milestone. Please adjust the "\
            "upper milestone(s) instead."
    
    elif modify_index == len(me.milestones) - 1: # modifying the highest
        # milestone
        lower_neighbor = me.milestones[modify_index].neighbors[0]
        assert radius > lower_neighbor.radius, "Cannot choose a radius "\
            "smaller than the lower adjacent milestone. Please adjust the "\
            "lower milestone(s) instead."
    else: # modifying an intermediate milestone
        upper_neighbor = me.milestones[modify_index].neighbors[1]
        lower_neighbor = me.milestones[modify_index].neighbors[0]
        assert radius > lower_neighbor.radius, "Cannot choose a radius "\
            "smaller than the lower adjacent milestone. Please adjust the "\
            "lower milestone(s) instead."
        assert radius < upper_neighbor.radius, "Cannot choose a radius "\
            "larger than the upper adjacent milestone. Please adjust the "\
            "upper milestone(s) instead."
            
    me.milestones[modify_index].radius = modify_radius
    return

def get_umbrella_info(me, milestone):
    absolute_directory = os.path.join(me.project.rootdir, milestone.directory)
    umbrella_glob = os.path.join(absolute_directory, 'md', 'umbrella', 
                                 'umbrella*.dcd')
    umbrella_dcds = glob.glob(umbrella_glob)
    try:
        catdcd_output = subprocess.check_output(['catdcd'] + umbrella_dcds)
    except FileNotFoundError:
        print('Make sure you have the program catdcd installed')
    except subprocess.CalledProcessError:
        catdcd_output = ''
    catdcd_lines = catdcd_output.splitlines()
    num_frames = 0
    for line in catdcd_lines:
        if line.startswith(b'Total'):
            num_frames = int(line.split()[-1])
    
    return num_frames
    
def get_fwd_rev_info(me, milestone):
    absolute_directory = os.path.join(me.project.rootdir, milestone.directory)
    reversal_glob = os.path.join(absolute_directory, 'md', 'fwd_rev', 
                                 'reverse*.dcd')
    reversal_dcds = glob.glob(reversal_glob)
    num_reversals = len(reversal_dcds)
    forward_glob = os.path.join(absolute_directory, 'md', 'fwd_rev', 
                                 'forward*.dcd')
    forward_dcds = glob.glob(forward_glob)
    num_forwards = len(forward_dcds)
    return num_reversals, num_forwards
    
def get_transition_info(me, milestone):
    absolute_directory = os.path.join(me.project.rootdir, milestone.directory)
    transition_filename = os.path.join(absolute_directory, 'md', 'fwd_rev', 
                                 'transition_info.xml')
    if not os.path.exists(transition_filename):
        return None
    else:
        trans_info = deserialize_transition_info(transition_filename)
        transition_strings = []
        incubation_time = trans_info['avg_incubation_time']
        for transition in trans_info['transitions']:
            source = transition['source']
            dest = transition['destination']
            count = transition['count']
            transition_strings.append('%d->%d: %d' % (source, dest, count))
        transition_string = ', '.join(transition_strings) \
            + ' incubation time: %.3f ps' % incubation_time
        return transition_string

def report_milestones(me):
    print("Milestones:")
    print("Index\tRadius\tNeighbors\tDirectory")
    for i, milestone in enumerate(me.milestones):
        if i == 0:
            print("%d\t%.2f\t" % (milestone.index, milestone.radius) , 
              milestone.neighbors[0].index, 
              '\t', milestone.directory)
        elif i == len(me.milestones)-1:
            print("%d\t%.2f\t" % (milestone.index, milestone.radius) , 
              milestone.neighbors[0].index,
              '\t', milestone.directory)
        else:
            print("%d\t%.2f\t" % (milestone.index, milestone.radius) , 
              milestone.neighbors[0].index, milestone.neighbors[1].index, 
              '\t', milestone.directory)
        umbrella_num_frames = get_umbrella_info(me, milestone)
        print('  Umbrella stage: %d frames found' % umbrella_num_frames)
        if umbrella_num_frames > 0:
            avg_distance = analyze.get_umbrella_avg_distance(me, milestone)
            print('  Umbrella average distance: %.2f angstroms' % avg_distance)
        fwd_rev_info = get_fwd_rev_info(me, milestone)
        print('  Fwd_rev Stage: %d reversal, %d forward dcd files found' \
        % fwd_rev_info)
        if not quick_analysis and (fwd_rev_info[0] > 0 or fwd_rev_info[1] > 0):
            fwd_rev_distances = analyze.get_fwd_rev_avg_distance(
                me, milestone)
            print('  Fwd_rev Stage: upward_distance: %.2f, ' \
                  'downward_distance: %.2f' % fwd_rev_distances)
        transition_info = get_transition_info(me, milestone)
        if transition_info is None:
            print('  Transitions: no data found')
        else:
            print('  Transitions: ' + transition_info)
        
        print('')

if __name__ == "__main__":
    print("Parse arguments")
    command = None
    directory = None
    
    assert len(sys.argv) >= 2, "At least one arguments required: 'command', "\
            "possibly 'seekrFile', and possibly {arguments...}."
    
    command = sys.argv[1]
    args = {}
    
    assert len(sys.argv) >= 3, "At least two arguments required: 'command', "\
        "and 'seekrFile', and, possibly {arguments...}."
    directory = sys.argv[2]
    
    assert command in ['add', 'modify', 'delete', 'report'], "Unknown command:"\
            " %s. Allowed commands: 'add', 'modify', 'delete', 'report'."
    if command == 'add':
        assert len(sys.argv) == 4, "Exactly one extra argument allowed for "\
                "'add' command: radius of adding milestone."
        args['radius'] = float(sys.argv[3])
    elif command == 'delete':
        assert len(sys.argv) == 4, "Exactly one extra argument allowed for "\
                "'delete' command: index of milestone to delete."
        args['index'] = int(sys.argv[3])
    elif command == 'modify':
        #pass # more complicated situation...
        assert len(sys.argv) == 5, "Exactly two extra arguments allowed for "\
                "'modify' command: index and radius of modified milestone."
        args['index'] = int(sys.argv[3])
        args['radius'] = float(sys.argv[4])
    elif command == 'report': # print a report about milestones in this system
        assert len(sys.argv) == 3, "No extra arguments allowed for 'report' "\
                "command."
    
    picklename = directory
    assert os.path.exists(picklename), "SEEKR pickle not found in provided "\
            "directory. Are you sure this is a SEEKR calculation directory?"
    
    print("Loading SEEKR calculation.")
    me = seekr.openSeekrCalc(picklename)
    
    # TODO: add support for multiple sites...
    
    if command == 'add':
        new_milestone, insert_index = add_milestone(me, args['radius'])
        # Save pickle and rename directories
        print("saving pickle...")
        add_milestone_directory(me, new_milestone, insert_index)
        me.save()
        print("Added milestone of radius", args['radius'], 
              "added to location:", insert_index)
    
    elif command == 'delete':
        del_dir, delete_index = delete_milestone(me, args['index'])
        delete_milestone_directory(del_dir, delete_index)
        me.save()
        print("Deleted milestone with index:", delete_index)
    
    elif command == 'modify':
        modify_milestone(me, args['index'], args['radius'])
        me.save()
        print("modified milestone", args['index'], "to have radius:", 
              args['radius'])
    
    elif command == 'report':
        report_milestones(me)
        
        
