'''
Created on May 8, 2018

@author: lvotapka

Main module for SEEKR runs
'''

from base import *
from milestones import *
from filetree import *
from building import *
import amber
import bd
from min_equil import *
from umbrella import *
from fwd_rev import *

def find_heavy_atoms(pdb_obj, index_range):
    '''Finds the heavy atoms within a subset of a PDB object.
    Input:
     - pdb_obj: a pdb2.py structure object that contains the atoms to search through.
     - index_range: a list of integers to search through.
    Output:
     - heavy_atoms: a list of integers that correspond to the indices of atoms
       that are heavy (non H).'''
    heavy_atoms = []
    for index in index_range:
        atom = pdb_obj.atoms[index]
        #print "atom.index:", atom.index, "index:", index, "element", atom.element
        if atom.element != 'H':
            heavy_atoms.append(index)
    return heavy_atoms
