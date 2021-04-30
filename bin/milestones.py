"""
Created on May 9, 2018

@author: lvotapka

Procedures for creating milestone objects
"""

import numpy as np
from math import sqrt
import cmath
import unittest


verbose = True

class Milestone_System():
  """An object for system information stored for each milestone."""
  def __init__(self):
    self.wet_holo_pdb_filename = ""
    self.dry_holo_pdb_filename = ""
    self.prmtop_filename = ""
    self.inpcrd_filename = ""
    self.umbrella_pdb_filename = ""
    self.system = None
    self.simulation = None

class Milestone():
  """Milestone superclass. It represents a surface in phase space that 
  is monitored for crossings.
  """
  def __init__(self, anchor, dimensions, index, siteid):
    self.fullname = ""
    self.directory = ""
    # the location of this milestone's anchor
    self.anchor = anchor
    #self.rotation = [1,0,0,0]
    self.neighbors = []
    # the index of the milestone
    self.index = index 
    # which set of milestones this one belongs to
    self.siteid = siteid 
    #self.center_type = center_type
    self.dimensions = dimensions
    # whether the milestone is kept stationary in space
    self.absolute = False
    # whether md simulations are run from this milestone
    self.md = False 
    # whether bd simulations are run from this milestone
    self.bd = False 
    # adjacent to a BD milestone
    self.bd_adjacent = None
    # a sink state milestone 
    self.end = False 
    self.openmm = Milestone_System()
    self.box_vectors = None
    
class Concentric_Spherical_Milestone(Milestone):
  """Concentric spherical milestones centered on an atom selection."""
  def __init__(self, index, siteid, absolute="False", md=True, bd=False):
    self.fullname = ""
    self.directory = ""
    # the location where the ligand was started
    self.anchor = None
    self.neighbors = []
    self.index = index
    self.siteid = siteid
    self.absolute = absolute
    self.md = md
    self.bd = bd
    self.bd_adjacent = None
    self.end = False
    # dimensions
    # the indices of the atoms in the system that define the 
    # center of this spherical milestone
    self.center_atom_indices = []
    # the x,y,z location of the center of this spherical milestone
    self.center_vec = None 
    self.radius = 0.0
    self.wet_holo_filename = ""
    self.dry_holo_filename = ""
    self.atom_selection_1 = None
    self.atom_selection_2 = None
    self.openmm = Milestone_System()
    self.config = None
    self.box_vectors = None
    
def quadratic_formula(a,b,c):
  """calculates the roots of the equation:
       a*x**2 + b*x + c = 0
     using the quadratic formula
  """
  results = []
  square_root = b*b - 4.0 * a * c
  for i in [-1.0, 1.0]:
    # then the roots will be real
    if square_root >= 0.0:
      root = (-b + (i * sqrt(square_root))) / (2.0 * a)
    # then the roots will be complex
    else: 
      root = (-b + (i * cmath.sqrt(square_root))) / (2.0 * a)
    results.append(root)
  return results

def find_anchor_from_vectors(origin, radius, vectors):
  """Given an origin and radius of the milestone, will follow the 
  vectors list by summing them until they cross the radius. Then, the 
  ligand will be placed at that location - the anchor.
  Input:
   - origin: 3-membered numpy array representing the center of the 
       spherical milestone
   - radius: float representing the radius of the milestone
   - vectors: a list of 3-membered arrays representing vectors that 
       define the
   'pathway' out of the receptor's active site
  Output:
   - anchor: 3-array representing the location of the anchor.
  """
  # starting from the very beginning of the
  summed_vectors = np.array([0.0,0.0,0.0])  
  for i in range(len(vectors)):
    cur_vector = vectors[i]
    cur_vector_prime = cur_vector / np.linalg.norm(cur_vector)
    # then this is the last vector, so the anchor has to be placed in 
    # line with this one
    if i == len(vectors)-1 or np.linalg.norm(summed_vectors + cur_vector) \
        > radius:
      # 'a' for the quadratic formula
      a = np.dot(cur_vector_prime,cur_vector_prime)
      # 'b' for the quadratic formula
      b = 2.0 * np.dot(summed_vectors,cur_vector_prime)
      # 'c' for the quadratic formula
      c = np.dot(summed_vectors,summed_vectors) - radius**2
      roots = quadratic_formula(a,b,c)
      # get the root with the smallest absolute value from the 
      # quadratic formula
      anchor = origin + summed_vectors + cur_vector_prime * min(map(abs,roots))
      return anchor
    summed_vectors += cur_vector
  return anchor
    
def generate_spherical_milestones(seekrcalc, atom_indices, origin, radius_list,
                                  siteid, vectors, absolute=False):
  """Create a list of concentric spherical milestones around an 
  atom selection.
  Input:
   - atom_indices: a list of integers representing the atom indices 
       that this milestone is centered on
   - origin: 3-array that should be equal (or approximately equal) to 
       the coordinates of the center of mass of the atom indices in 
       the input structure
   - radius_list: a list of floats representing the radii of the 
       concentric spheres in the list of milestones
   - siteid: integer representing the index of this list of milestones
   - vectors: a list of lists of floats, or of numpy arrays, that 
       define the pathway out of the active site, and in turn, the 
       locations where the ligand will be placed
   - absolute: boolean defining whether this milestone will remain 
       stationary in space regardless of how the receptor moves  
  Output:
   - milestones: a list of milestone objects
  """
    
  milestones = []
  # the lowest radius
  lowest_radius = radius_list[0]
  #TODO: remove
  #total_spherical_milestones = 0
  spheres_in_site = len(radius_list)
  # TODO: remove
  #cur_vector = vectors[0]
  # loop through all the milestone radii
  for i in range(spheres_in_site):
    radius = radius_list[i]
    # now we need to find the anchor - the location where 
    # the ligand will be placed
    anchor = find_anchor_from_vectors(origin, radius, vectors)
    anchor_tuple = tuple(anchor)
    milestone = Concentric_Spherical_Milestone(i, siteid, absolute, 
                                               md=seekrcalc.project.md, 
                                               bd=seekrcalc.project.bd)
    milestone.center_atom_indices = atom_indices
    milestone.center_vec = origin
    milestone.radius = radius
    milestone.anchor = anchor
    if i > 0:
      # if we are the furthest milestone in either direction, 
      # then include no neighbors
      milestone.neighbors.append(i-1)
    #else:
      #if not k_off: milestone.end = True
    if i < spheres_in_site - 1:
      milestone.neighbors.append(i+1)
      milestone.bd = False
    else: # then this is the outermost one, so set BD to true
      milestone.bd = True
      # make the adjacent milestone the previously created md milestone
      milestone.bd_adjacent = milestones[-1]
      milestone.end = True
      milestone.md = False
    milestone.fullname = "%d_%d_%s_%.1f_%.1f_%.1f" \
        % (i+(spheres_in_site*int(siteid)), i, siteid, anchor[0], 
           anchor[1], anchor[2])
    milestones.append(milestone)
  return milestones

def print_spherical_milestone_info(milestones):
  """Print the information in all milestones for debug and 
  verbosity purposes.
  Input:
   - milestones: a list of milestone objects to print information about.
  Output:
   - None
   """
  for milestone in milestones:
    print "Milestone index:", milestone.index, "siteid:", milestone.siteid, \
        "origin:", milestone.center_vec, \
        "center_atoms:", milestone.center_atom_indices, 
    print "  num neighbors:", len(milestone.neighbors), \
        "anchor:", milestone.anchor, "md:", milestone.md, "bd:", milestone.bd,\
        "bd_adjacent:", milestone.bd_adjacent
    print "  end:", milestone.end, "fullname:", milestone.fullname, \
        "radius:", milestone.radius
  return

class Test_milestones(unittest.TestCase):
  # several test cases to ensure the functions in this module are 
  # working properly
  # test this function
  def test_main(self): 
    print "WARNING: this module does not have comprehensive unit tests"
    
  def test_quadratic_formula(self):
    # test negative roots
    self.assertEqual(quadratic_formula(1.0,4.0,3.0), [-3.0,-1.0]) 
    # test positive roots
    self.assertEqual(quadratic_formula(3.0,-13.0,12.0), [4.0/3.0,3.0]) 
     # test zero roots
    self.assertEqual(quadratic_formula(1.0,-2.0,0.0), [0.0,2.0])
    # test zero roots
    self.assertEqual(quadratic_formula(1.0,0.0,0.0), [0.0,0.0]) 
    # test imaginary roots
    self.assertEqual(quadratic_formula(1.0,0.0,4.0), 
                     [complex(0,-2.0),complex(0,2.0)]) 
    return

  def test_find_rad_of_skewed_vectors(self):
    A = np.array([0.0,1.0,0.0])
    B = np.array([1.0,0.0,0.0])
    B2 = np.array([1.0,1.0,0.0])
    origin = np.array([0,0,0])
    """ # TODO: complete this later
    # check a vector against itself
    self.assertEqual(find_anchor_from_vectors(origin, A, 2.0), 1.0) 
    # check a vector against right angle
    self.assertEqual(find_anchor_from_vectors(A, B, 2.0), sqrt(3.0)) 
    # check a vector against right angle
    self.assertEqual(find_anchor_from_vectors(2.0*A, B*2.5, 3.0), sqrt(5.0)) 
    # check a vector against 45 degree angle
    self.assertEqual(find_anchor_from_vectors(B2, A, 2.0), sqrt(3.0)-1.0) 
    """
    return
