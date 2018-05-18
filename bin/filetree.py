#!/usr/bin/python

'''
filetree.py

Contains functions to be called by seekr that generate the necessary filetree to execute simulations.

'''
import sys, os, math, shutil, subprocess #, make_fxd
import unittest
import numpy as np
import pdb2 as pdb

verbose = True
mdtree = {'building':{},'min':{},'temp_equil':{}, 'umbrella':{}, 'fwd_rev':{}} 
bdtree = {'ens':{}}


class Filetree():
  '''defines a file tree: a framework of directories to be populated with files.'''
  def __init__(self, tree):
    '''
    Input:
        tree: a nested dictionary object defining the file structure to create.
            example: {'file1':{}, 'file2':{'file3':{}}}
    Output: 
        None
    '''
    self.tree = tree
    return

  def make_tree(self, rootdir, branch={}):
    '''will construct the file tree given a root directory, populating
it with all subsequent branches and leaves in the file tree.
    Input: 
        rootdir: a string indicating the path in which to place the root of the file tree
        branch: optional nested dictionary structure of the tree. (see Filetree.__init__)
    Output:
        None
    '''
    assert os.path.isdir(rootdir), "rootdir argument must be a real directory"
    if not branch: branch = self.tree
    for subbranch in branch.keys():
      # first create each subbranch
      subbranchpath = os.path.join(rootdir,subbranch)
      if not os.path.exists(subbranchpath):
        os.mkdir(subbranchpath)
      if not branch[subbranch]: # if its an empty list, then we have a leaf
        continue
      else: # then we can descend further, recursively
        self.make_tree(subbranchpath, branch=branch[subbranch])
    return

def generate_filetree(seekrcalc):
  '''Generate the filetree for a SEEKR run within seekrcalc.project.rootdir directory.
      Input:
          seekrcalc: a Seekr_Calculation object of commonly-used settings in SEEKR runs.
      Output:
          None
  '''
  rootdir = seekrcalc.project.rootdir
  if seekrcalc.project.empty_rootdir and os.path.isdir(seekrcalc.project.rootdir): # then we're in test mode, delete anything that is in the rootdir
    print "'empty_rootdir' set to True. Deleting all subdirectories/files in rootdir:", seekrcalc.project.rootdir
    shutil.rmtree(seekrcalc.project.rootdir) # remove all the files within the grid subdirectory
  if not os.path.isdir(seekrcalc.project.rootdir): # if the rootdir doesn't exist yet
    os.mkdir(seekrcalc.project.rootdir) # then create it
  
  config_dirlist = []
  config_counter = 0
  for i, milestone in enumerate(seekrcalc.milestones): # for each configuration, it gets its own directory
    #wet_config = wet_configs[i]
    anchor_name = "anchor_%s" % (milestone.fullname)
    anchor_filetree = Filetree({anchor_name:{}})
    if verbose: print "Creating directory:", anchor_name
    anchor_filetree.make_tree(rootdir) # create this anchor's directory
    milestone.directory = anchor_name # update this milestones directory information
    anchor_dir = os.path.join(rootdir, anchor_name)
    
    # MD filetree
    md_file_path={}
    if seekrcalc.project.md == True and milestone.md == True: # then prep this anchor for an MD simulation
      md_dir=os.path.join(anchor_dir,'md') # directory for MD
      md_filetree=Filetree({'md':mdtree})
      md_filetree.make_tree(anchor_dir) # make the MD filetree
      for key in mdtree.keys():
        md_file_path[key] = os.path.join(md_dir,key)
      wet_holo_filename=os.path.join(md_dir,'holo_wet.pdb')
      #wet_config.save(wet_holo_filename, amber=True, standard=False) # write the holo structure into the md directory
      #md_file_path['wet_holo'] = wet_holo_filename
      milestone.openmm.wet_holo_pdb_filename = wet_holo_filename
    
    # BD filetree
    bd_file_path={}
    bd_dir = None
    if seekrcalc.project.bd == True and milestone.bd == True: # then prep this anchor for BD simulation
      #dry_config= wet_config #pdb.dry(wet_config) # since we didn't iterate through this list, we must use the index
      bd_dir=os.path.join(anchor_dir,'bd') # directory for MD
      bd_filetree=Filetree({'bd':bdtree})
      bd_filetree.make_tree(anchor_dir) # make the MD filetree
      for key in bdtree.keys():
        bd_file_path[key] = os.path.join(bd_dir,key)
      dry_holo_filename=os.path.join(bd_dir,'holo_dry.pdb')
      #dry_config.save(dry_holo_filename, amber=True, standard=False) # write the holo structure into the md directory
      #md_file_path['dry_holo'] = dry_holo_filename
      milestone.openmm.dry_holo_pdb_filename = dry_holo_filename

    config_dirlist.append(anchor_name)
    #if md_file_path: md_file_paths.append(md_file_path)
    #if bd_dir: bd_file_paths.append(bd_dir) # bd_file_path
    config_counter+=1 # increment the loop counter
  
  #seekrcalc.building.md_file_paths = md_file_paths
  #seekrcalc.building.bd_file_paths = bd_file_paths
  seekrcalc.building.config_dirlist = config_dirlist
  return 

class Test_filetree_functions(unittest.TestCase):
  # several test cases to ensure the functions in this module are working properly
  def test_main(self): # test this function
    pass # Skipping because this function requires complicated objects to run

  def test_Filetree(self): # test the Filetree class
    # test init
    empty_filetree = {'mytest':{}}
    test = Filetree(empty_filetree)
    self.assertEqual(test.tree, empty_filetree)
    # test make_tree
    test.make_tree('/tmp')
    test_filetree = {'test1':{}, 'test2':{'test3':{}}}
    test_rootdir = '/tmp/mytest'
    tree=Filetree(test_filetree)
    tree.make_tree(test_rootdir) # make the test filetree
    self.assertTrue(os.path.exists('/tmp/mytest/test1'))
    self.assertTrue(os.path.exists('/tmp/mytest/test2')) # test to see whether these directories have been properly created
    self.assertTrue(os.path.exists('/tmp/mytest/test2/test3'))


if __name__ == "__main__":
  print "Running unit tests for filetree.py"
  unittest.main() # run tests of all functions
