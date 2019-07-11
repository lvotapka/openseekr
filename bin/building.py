'''
Created on May 9, 2018

@author: lvotapka

Module for generating starting structures, topology files, and coordinate files
for later minimization, etc.
'''
import pdb2 as pdb # custom library for reading/writing pdb files
import numpy as np
from copy import deepcopy # needed to keep track of separate structure objects
import time
import os
import unittest
import cPickle as pickle

verbose = True # whether stuff is printed in detail

def structures_clash(structure1, structure2, tolerance=0.0):
  '''determines whether two structures have any clashing atoms
    the smallest structure should be the first argument to allow
    this function to run quickly as possible
    Input:
     - structure1: a Structure() object, probably representing the receptor
     - structure2: a Structure() object, probably representing the ligand
     - tolerance: float, The distance, in Angstroms, to subtract from the addition of 
         the radii before a clash is rejected.
    Output:
     - clashing: a boolean representing whether the structures are clashing or 
         not.
  '''
  if verbose: print "searching for clashes between structures", structure1.struct_id, "and", structure2.struct_id #If true execute the print statement.
  struct1_com = pdb.center_of_mass(structure1) # find center of mass of structure1
  struct1_rad = pdb.molecular_radius(structure1) # find molecular radius
  clashing = False

  for atom2 in structure2.atoms: # for every atom in structure2
    if atom2.radius == '0.0':
      atom2_radius = pdb.radii[atom2.element] # if not included, then retrieve a standard radius for this atom
      if verbose: print "using dictionary on atom2. atom2.radius:", atom2.radius
    else:
      atom2_radius = float(atom2.radius)

    # the statement below saves some computation time by being a semi-divide and conquer method. Could be better but I'm sure it will work just fine
    if np.linalg.norm(np.array(struct1_com) - np.array(atom2.coords)) < struct1_rad + atom2_radius:
      # then we have to check this atom against every struct1 atom
      for atom1 in structure1.get_atoms(): # for every atom in structure1
        if atom1.radius == '0.0':
          atom1_radius = pdb.radii[atom1.element] # if not included, then retrieve a standard radius for this atom
          if verbose: print "using dictionary on atom1. atom1.radius:", atom1.radius
        else:
          atom1_radius = float(atom1.radius)
        atom_dist = atom2_radius + atom1_radius
        if np.linalg.norm(np.array(atom1.coords) - np.array(atom2.coords)) < atom_dist - tolerance:
          return True
  return clashing

def pickle_or_load(filename, picklename, parser, struc_name="pickle",pqr=False):
  '''for large files, instead of parsing, they can be saved and loaded much more 
  quickly as a pickle. 
  Input:
   - filename: string, the name of the pqr file to save
   - picklename: string, the name of the pickle file to check, and perhaps load
   - parser: a Big_PDBParser() object used to parse the pqr file
   - struct_name: string, the name of the structure object
   - pqr: boolean, whether the file should be opened as a pqr file
  Output:
   - The Structure() object of the loaded molecule
  '''
  if os.path.exists(picklename) and os.path.getmtime(picklename) > os.path.getmtime(filename): # if the pickle has been most recently modified
    # load the pickle
    if verbose: print "reading pickle:", picklename
    our_file=open(picklename, 'rb')
    our_obj=pickle.load(our_file)
    our_file.close()
  else:
    # then load the file itself and save the pickle
    our_obj=parser.get_structure(struc_name, filename, pqr=pqr, conventional=False) # load the file
    if verbose: print "writing pickle:", picklename
    our_file=open(picklename, 'wb')
    pickle.dump(our_obj, our_file, protocol=-1) # protoco=-1 means to use the fastest available protocol method
    our_file.close()
  return our_obj

def generate_configs(seekrcalc):
  '''Generates the apo and holo configurations of the receptors and ligands.
  Input:
   - seekrcalc: a SeekrCalculation() object that contains all the settings for
       the SEEKR calculation
  Output:
   - None
  '''
  if verbose: print "Generating structural configurations..."
  parser = pdb.Big_PDBParser()
  if verbose: print "now loading structures"
  holo_config_wet = None
  # Read and/or create pickle files for the structures to save I/O time
  ligand_pkl_filename = os.path.join(seekrcalc.project.rootdir, "ligand.pkl")
  receptor_pkl_wet_filename = os.path.join(seekrcalc.project.rootdir, "receptor.pkl")
  #receptor_pkl_dry_filename = os.path.join(seekrcalc.project.rootdir, "receptor_dry.pkl")
  receptor_pkl_dry_pqr_filename = os.path.join(seekrcalc.project.rootdir, "receptor_dry_pqr.pkl")

  ligand=pickle_or_load(seekrcalc.building.lig_dry_pqr_filename, ligand_pkl_filename, parser, struc_name="ligand", pqr=True)
  receptor_wet=pickle_or_load(seekrcalc.building.rec_wet_pdb_filename, receptor_pkl_wet_filename, parser, struc_name="receptor_wet", pqr=False)
  receptor_dry=pickle_or_load(seekrcalc.building.rec_dry_pqr_filename, receptor_pkl_dry_pqr_filename, parser, struc_name="receptor_dry_pqr", pqr=True)
  
  seekrcalc.building.ligand = ligand
  seekrcalc.building.rec_wet_pdb_filename = receptor_wet
  seekrcalc.building.rec_dry_pqr_filename = receptor_dry
  milestones = seekrcalc.milestones
  configs = []
  starttime = time.time()
  if verbose: print "Generating", len(milestones), "ligand configurations"
  for milestone in milestones:
    new_ligand = deepcopy(ligand) # copy the entire structure
    if milestone.lig_center is not None:
      struct_center = milestone.lig_center
    else:
      struct_center = pdb.center_of_mass(ligand)
    new_ligand.moveby(-struct_center) # reposition ligand structure over the origin
    new_ligand.moveby(milestone.anchor) # move ligand to the anchor location
    new_ligand.struct_id = milestone.fullname
    configs.append(new_ligand)
  endtime = time.time()
  if verbose:
    "Generate_configs complete. Total number: %d. Elapsed time: %d s" % (len(configs), endtime-starttime,)
    print "Now running through all configurations to combine ligand with receptor."
    if seekrcalc.building.reject_clashes:
      print "Rejecting steric clashes..."
    else:
      print "Ignoring steric clashes..."
  
  if seekrcalc.project.bd:
    seekrcalc.browndye.starting_lig_config = configs[0]
  
  for i, milestone in enumerate(milestones): # run thru all the configurations of ligands
    lig_config = configs[i]
    milestone.config = deepcopy(lig_config)
    if seekrcalc.building.reject_clashes:
      if structures_clash(lig_config, receptor_dry, tolerance=.00000000001): continue # if there's a clash
    if milestone.md:
      #pdb.TER_resnames.append(seekrcalc.building.lig_resname) # TODO: marked for removal
      holo_config_wet, insert_index, last_ligand_index = pdb.ligmerge(lig_config, receptor_wet, verbose=False)
      holo_config_wet.struct_id = lig_config.struct_id # set the structure description to the same as the ligand
      holo_config_wet.renumber_indeces() # to number the indeces consecutively
      wet_holo_filename = milestone.openmm.wet_holo_pdb_filename
      if verbose: print "writing file:", wet_holo_filename
      holo_config_wet.save(wet_holo_filename, amber=True, standard=False) # write the holo structure into the md directory
      last_insert_index = insert_index
      last_last_ligand_index = last_ligand_index
    if milestone.bd and seekrcalc.project.bd:
      holo_config_dry, insert_index, last_ligand_index = pdb.ligmerge(lig_config, receptor_dry, verbose=False)
      holo_config_dry.struct_id = lig_config.struct_id # set the structure description to the same as the ligand
      holo_config_dry.renumber_indeces() # to number the indeces consecutively
      dry_holo_filename = milestone.openmm.dry_holo_pdb_filename
      if verbose: print "writing file:", dry_holo_filename
      holo_config_dry.save(dry_holo_filename, amber=True, standard=False) # write the holo structure into the md directory
  
  assert holo_config_wet != None, "All structures clashed. Consider adding more milestones. Exiting..."
  return holo_config_wet, last_insert_index, last_last_ligand_index

