'''
Created on Mar 15, 2019

@author: lvotapka

Module for modifying, adding, or deleting milestones
'''

import seekr
#from seekr import amber
import sys, os, math, shutil
import numpy as np

print "Parse arguments"
command = None
directory = None

assert len(sys.argv) >= 3, "At least two arguments required: 'command', and 'directory', and, possibly {arguments...}."

command = sys.argv[1]
directory = sys.argv[2]
args = {}

assert command in ['add', 'modify', 'delete', 'report'], "Unknown command: %s. Allowed commands: 'add', 'modify', 'delete', and 'report'."
if command == 'add':
  assert len(sys.argv) == 4, "Only one extra argument allowed for 'add' command: radius of adding milestone."
  args['radius'] = float(sys.argv[3])
elif command == 'delete':
  assert len(sys.argv) == 4, "Only one extra argument allowed for 'delete' command: index of milestone to delete."
  args['index'] = int(sys.argv[3])
elif command == 'modify':
  #pass # more complicated situation...
  assert len(sys.argv) == 5, "Exactly two extra arguments allowed for 'modify' command: index and radius of modified milestone."
  args['index'] = int(sys.argv[3])
  args['radius'] = float(sys.argv[4])
elif command == 'report': # print a concise report about all of the milestones in this system
  assert len(sys.argv) == 3, "No extra arguments allowed for 'report' command."
  
  
picklename = os.path.join(directory, 'seekr_calc.pickle')
assert os.path.exists(picklename), "SEEKR pickle not found in provided directory. Are you sure this is a SEEKR calculation directory?"

print "Loading SEEKR calculation."
me = seekr.openSeekrCalc(picklename)

# TODO: add support for multiple sites...
siteid = 0
spheres_in_site = 0

if command == 'add':
  insert_index = None # index to insert before
  for i, milestone in enumerate(me.milestones):
    if milestone.radius > args['radius'] and insert_index == None:
      insert_index = i
      absolute = milestone.absolute
      atom_indices = milestone.center_atom_indices
      center_vec = milestone.center_vec
      
    if milestone.siteid == siteid:
      spheres_in_site += 1
  
  assert insert_index != None, 'The provided radius exists beyond the current BD milestone. This is not currently supported. Move the existing BD milestone, and add a new milestone below it.'
  
  # create the new milestone
  new_milestone = seekr.Concentric_Spherical_Milestone(i, siteid, absolute, md=me.project.md, bd=me.project.bd)
  anchor = np.array([0.0, 0.0, 0.0]) # TODO: may need to be modified in the future
  #anchor = seekr.find_anchor_from_vectors(center_vec, args['radius'], vectors)
  new_milestone.index = insert_index
  new_milestone.center_atom_indices = atom_indices
  new_milestone.center_vec = center_vec
  new_milestone.radius = args['radius']
  new_milestone.anchor = anchor
  lower_neighbor_index = None
  upper_neighbor_index = None
  if insert_index > 0:
    lower_neighbor_index = insert_index - 1
    new_milestone.neighbors.append(lower_neighbor_index)
  
  upper_neighbor_index = insert_index + 1
  new_milestone.neighbors.append(upper_neighbor_index)
  new_milestone.bd = False
  new_milestone.fullname = "%d_%d_%d_%.1f_%.1f_%.1f" % (insert_index+(spheres_in_site*int(siteid)), insert_index, siteid, anchor[0], anchor[1], anchor[2])
  new_milestone.directory = "anchor_%s" % new_milestone.fullname
  
  # modify neighbor milestones
  # Lower neighbor doesn't need modification
  for i in range(insert_index, len(me.milestones)): # all upper milestones need modifying
    me.milestones[i].index = i+1
    me.milestones[i].neighbors[0] = i
    if len(me.milestones[i].neighbors) == 2:
      me.milestones[i].neighbors[1] = i+2
    # parse this milestone's name to construct a new one
    old_fullname = me.milestones[i].fullname
    old_dirname = me.milestones[i].directory
    id1, id2, mysiteid, anchor0, anchor1, anchor2 = me.milestones[i].fullname.split('_')
    new_fullname = "%d_%d_%s_%s_%s_%s" % (int(id1)+1, int(id2)+1, mysiteid, anchor0, anchor1, anchor2)
    new_dirname = 'anchor_%s' % new_fullname
    me.milestones[i].fullname = new_fullname
    me.milestones[i].directory = new_dirname
    old_path = os.path.join(me.project.rootdir, old_dirname)
    old_path_tmp = old_path+"temp"
    new_path = os.path.join(me.project.rootdir, new_dirname)
    new_path_tmp = new_path+"temp"
    if me.milestones[i].openmm.prmtop_filename:
      me.milestones[i].openmm.prmtop_filename = os.path.join(new_path, 'md', 'building','holo.parm7')
    if me.milestones[i].openmm.inpcrd_filename:
      me.milestones[i].openmm.inpcrd_filename = os.path.join(new_path, 'md', 'building','holo.rst7')
    if os.path.exists(new_path):
      #tmp_name = os.path.join(me.project.rootdir, 'temp_dir')
      if os.path.exists(old_path_tmp):
        print "temporarily renaming", old_path_tmp, "to", new_path_tmp
        os.rename(old_path_tmp, new_path_tmp)
      else:
        print "temporarily renaming", old_path, "to", new_path_tmp
        os.rename(old_path, new_path_tmp)
    else:
      if os.path.exists(old_path_tmp):
        print "renaming files:", old_path_tmp, "to", new_path
        os.rename(old_path_tmp, new_path)
      else:
        print "renaming files:", old_path, "to", new_path
        os.rename(old_path, new_path)

  
  # insert new milestone into the milestone list in order
  me.milestones.insert(insert_index, new_milestone)
  anchor_dir = os.path.join(me.project.rootdir, new_milestone.directory)
  anchor_filetree = seekr.Filetree({new_milestone.directory:{}})
  print "Making directory:", anchor_dir
  anchor_filetree.make_tree(me.project.rootdir) # create this anchor's directory
  md_dir=os.path.join(anchor_dir,'md') # directory for MD
  md_filetree=seekr.Filetree({'md':seekr.mdtree})
  md_filetree.make_tree(anchor_dir) # make the MD filetree
  print "Making directory:", md_dir
  
  
  # Save pickle and rename directories
  print "saving pickle..."
  me.save()
  print "New milestone of radius", args['radius'], "added."
  

elif command == 'delete':
  delete_index = args['index']
  for i in range(delete_index+1, len(me.milestones)): # all upper milestones need modifying
    me.milestones[i].index = i-1
    me.milestones[i].neighbors[0] = i-2
    if len(me.milestones[i].neighbors) == 2:
      me.milestones[i].neighbors[1] = i
    # parse this milestone's name to construct a new one
    old_fullname = me.milestones[i].fullname
    old_dirname = me.milestones[i].directory
    id1, id2, mysiteid, anchor0, anchor1, anchor2 = me.milestones[i].fullname.split('_')
    new_fullname = "%d_%d_%s_%s_%s_%s" % (int(id1)-1, int(id2)-1, mysiteid, anchor0, anchor1, anchor2)
    new_dirname = 'anchor_%s' % new_fullname
    #print "old_dirname:", old_dirname, "new_dirname:", new_dirname
    me.milestones[i].fullname = new_fullname
    me.milestones[i].directory = new_dirname
    old_path = os.path.join(me.project.rootdir, old_dirname)
    new_path = os.path.join(me.project.rootdir, new_dirname)
    print "renaming", old_path, "to", new_path
    os.rename(old_path, new_path)
    
    
  del_milestone = me.milestones.pop(delete_index)
  print "deleting milestone index:", del_milestone.index, "radius:", del_milestone.radius
  del_dir = os.path.join(me.project.rootdir, del_milestone.directory)
  shutil.rmtree(del_dir)
  
  me.save()

elif command == 'modify':
  modify_index = args['index']
  modify_radius = args['radius']
  me.milestones[modify_index].radius = modify_radius
  print "modifying milestone", modify_index, "to have radius:", modify_radius
  me.save()

elif command == 'report':
  print "Milestones:"
  print "Index\tRadius\tNeighbors\tDirectory"
  for i, milestone in enumerate(me.milestones):
    print "%d\t%f\t" % (milestone.index, milestone.radius) , milestone.neighbors, '\t', milestone.directory
    
