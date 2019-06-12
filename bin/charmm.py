'''
Merge ligand and receptor files to write parameter/topology/coordinate files
for Charmm simulations.

Created on Jun 11, 2019

@author: lvotapka
'''

import os, string
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from copy import deepcopy
import parmed
import psfmerge

verbose = True

class CharmmSettings():
  def __init__(self):
    self.receptor_psf = ''
    self.ligand_psf = ''
    self.insert_index = -1

def charmm_building(seekrcalc, milestone, charmm_settings):
  '''the pre-minimization procedure for charmm ff simulations
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to prepare the AMBER system for
   ###- amber_settings: the AmberSettings() object that defines settings for the
       ###tleap calculation.
  Output:
   - None
  '''
  i = milestone.index
  building = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'building')
  pdb_filename = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'holo_wet.pdb')
  psf_filename = os.path.join(building, 'holo_wet.psf')
  
  if not os.path.exists(pdb_filename):
    print "Skipping milestone %d: no holo structure detected at %s." % (i, pdb_filename)
    return
  
  insert_index = charmm_settings.insert_index
  newpsf = psfmerge.merge_psf_files(charmm_settings.receptor_psf, 
                                    charmm_settings.ligand_psf, 
                                    insert_index = insert_index) # merge the psf files
  newpsf.write(psf_filename) # write the new psf files
  return


  
  
