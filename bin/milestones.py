'''
Created on June 11, 2019


@author: lvotapka
@author: astokely

Procedures for creating spherical and planar milestone objects
'''

import numpy as np
from math import sqrt
import cmath
import unittest
from copy import deepcopy


verbose = True

class Milestone_System():
  '''An object for system information stored for each milestone.'''
  def __init__(self):
    self.wet_holo_pdb_filename = ''
    self.dry_holo_pdb_filename = ''
    self.prmtop_filename = '' # AMBER inputs
    self.inpcrd_filename = ''
    self.psf_filename = ''
    #self.rtf_filename = '' # CHARMM inputs # TODO: remove
    self.par_filename = ''
    self.charmm_params_filename_list = []
    self.umbrella_pdb_filename = ''
    self.system = None
    self.simulation = None

class Milestone():
  '''Milestone superclass. It represents a surface in phase space that is monitored for crossings'''
  def __init__(self, anchor, dimensions, index, siteid):
    self.type = ""
    self.fullname = ""
    self.directory = ""
    self.anchor = anchor # the location of this milestone's anchor
    #self.rotation = [1,0,0,0]
    self.neighbors = []
    self.index = index # the index of the milestone
    self.siteid = siteid # which set of milestones this one belongs to
    #self.center_type = center_type
    self.dimensions = dimensions
    self.absolute = False # whether the milestone is kept stationary in space
    self.md = False # whether md simulations are run from this milestone
    self.bd = False # whether bd simulations are run from this milestone
    self.bd_adjacent = None # adjacent to a BD milestone
    self.end = False # a sink state milestone
    self.openmm = Milestone_System()
    self.box_vectors = None
    
class Concentric_Spherical_Milestone(Milestone):
  '''Concentric spherical milestones centered on an atom selection.'''
  def __init__(self, index, siteid, absolute='False', md=True, bd=False):
    self.type = 'spherical'
    self.fullname = ''
    self.directory = ''
    self.anchor = None # the location where the ligand was started
    self.neighbors = []
    self.index = index
    self.siteid = siteid
    self.absolute = absolute
    self.md = md
    self.bd = bd
    self.bd_adjacent = None
    self.end = False
    # dimensions
    self.center_atom_indices = [] # the indices of the atoms in the system that define the center of this spherical milestone
    self.center_vec = None # the x,y,z location of the center of this spherical milestone
    self.lig_center = None
    self.radius = 0.0
    self.wet_holo_filename = ''
    self.dry_holo_filename = ''
    self.atom_selection_1 = None # TODO: review role of this selection vs. center_atom_indices
    self.atom_selection_2 = None
    self.openmm = Milestone_System()
    self.config = None
    self.box_vectors = None
    
class Planar_Z_Milestone(Milestone):
  '''Planar milestones perpendicular to Z-axis and centered on an atom selection.'''
  def __init__(self, index, siteid, absolute='False', md=True, bd=False):
    self.type = 'planar_z'
    self.fullname = ''
    self.directory = ''
    self.anchor = None # the location where the ligand was started
    self.neighbors = []
    self.index = index
    self.siteid = siteid
    self.absolute = absolute
    self.md = md
    self.bd = bd
    self.bd_adjacent = None
    self.end = False
    # dimensions
    self.center_atom_indices = [] # the indices of the atoms in the system that define the center of this spherical milestone
    self.center_vec = None # the x,y,z location of the center of this spherical milestone
    self.lig_center = None
    self.offset = 0.0
    self.wet_holo_filename = ''
    self.atom_selection_1 = None
    self.atom_selection_2 = None
    self.openmm = Milestone_System()
    self.config = None
    self.box_vectors = None
    
class Ellipsoidal_Milestone(Milestone):
  '''Ellipsoidal milestones centered on an atom selection.'''
  def __init__(self, index, siteid, absolute='False', md=True, bd=False):
    self.type = 'ellipsoidal'
    self.fullname = ''
    self.directory = ''
    self.anchor = None # the location where the ligand was started
    self.neighbors = []
    self.index = index
    self.siteid = siteid
    self.absolute = absolute
    self.md = md
    self.bd = bd
    self.bd_adjacent = None
    self.end = False
    # dimensions
    self.locus_1_atom_indices = [] # the indices of the atoms in the system that define the center of this spherical milestone
    self.locus_1_vec = None # the x,y,z location of the center of this spherical milestone
    self.locus_2_atom_indices = [] # the indices of the atoms in the system that define the center of this spherical milestone
    self.locus_2_vec = None # the x,y,z location of the center of this spherical milestone
    self.lig_center = None
    self.nu = 0.0
    self.wet_holo_filename = ''
    self.dry_holo_filename = ''
    self.atom_selection_locus_1 = None
    self.atom_selection_locus_2 = None
    self.atom_selection_ligand = None
    self.openmm = Milestone_System()
    self.config = None
    self.box_vectors = None
    
class RMSD_Milestone(Milestone):
  '''RMSD milestones centered on a list of atom selections.'''
  def __init__(self, index, siteid, absolute='False', md=True, bd=False):
    self.type = 'rmsd'
    self.fullname = ''
    self.directory = ''
    self.anchor = None # the location where the ligand was started
    self.neighbors = []
    self.index = index
    self.siteid = siteid
    self.absolute = absolute
    self.md = md
    self.bd = bd
    self.bd_adjacent = None
    self.end = False
    # dimensions
    self.atom_indices1 = [] # the indices of the atoms in the system that define the center of this spherical milestone
    self.atom_points1 = []
    self.atom_vec1 = None # the x,y,z location of the center of this spherical milestone
    self.atom_indices2 = []
    self.atom_points2 = []
    self.atom_vec2 = None
    self.lig_center = None
    self.radius = 0.0
    self.wet_holo_filename = ''
    self.dry_holo_filename = ''
    self.openmm = Milestone_System()
    self.config = None
    self.box_vectors = None
    
def quadratic_formula(a,b,c):
  '''calculates the roots of the equation:
       a*x**2 + b*x + c = 0
     using the quadratic formula
  '''
  results = []
  square_root = b*b - 4.0 * a * c
  for i in [-1.0, 1.0]:
    if square_root >= 0.0: # then the roots will be real
      root = (-b + (i * sqrt(square_root))) / (2.0 * a)
    else: # then the roots will be complex
      root = (-b + (i * cmath.sqrt(square_root))) / (2.0 * a)
    results.append(root)
  return results

def find_spherical_anchor_from_vectors(origin, radius, vectors):
  '''Given an origin and radius of the milestone, will follow the vectors list
  by summing them until they cross the radius. Then, the ligand will be placed
  at that location - the anchor.
  Input:
   - origin: 3-membered numpy array representing the center of the spherical 
  milestone
   - radius: float representing the radius of the milestone
   - vectors: a list of 3-membered arrays representing vectors that define the
   'pathway' out of the receptor's active site
  Output:
   - anchor: 3-array representing the location of the anchor.
  '''
  summed_vectors = np.array([0.0,0.0,0.0]) # starting from the very beginning of the origin 
  for i in range(len(vectors)):
    cur_vector = vectors[i]
    cur_vector_prime = cur_vector / np.linalg.norm(cur_vector)
    if i == len(vectors)-1 or np.linalg.norm(summed_vectors + cur_vector) > radius: # then this is the last vector, so the anchor has to be placed in line with this one
      a = np.dot(cur_vector_prime,cur_vector_prime) # 'a' for the quadratic formula
      b = 2.0 * np.dot(summed_vectors,cur_vector_prime) # 'b' for the quadratic formula
      c = np.dot(summed_vectors,summed_vectors) - radius**2 # 'c' for the quadratic formula
      roots = quadratic_formula(a,b,c)
      anchor = origin + summed_vectors + cur_vector_prime * min(map(abs,roots)) # get the root with the smallest absolute value from the quadratic formula
      return anchor
    summed_vectors += cur_vector
  return anchor

def get_points_diff(points1, points2):
  '''Find vectors (points_diff) that points from points1 to points2'''
  assert len(points1) ==  len(points2), "the lengths of points1 and points2 must be equal."
  n = len(points1)
  points_diff = []
  for i in range(n):
    points_diff.append(points1[i] - points2[i])
  return points_diff

def points_rmsd(points_diff):
  '''Finds the RMSD between two sets of points'''
  msd = 0.0
  n = len(points_diff)
  for i in range(n):
    #point1 = points1[i]
    #point2 = points2[i]
    #points_diff = point2 - point1
    msd += np.dot(points_diff[i], points_diff[i])
  rmsd = np.sqrt(msd)
  return rmsd

def get_points_com(points):
  '''compute the center of mass of a vector.'''
  com = np.array([0.0, 0.0, 0.0])
  n = len(points)
  for i in range(n):
    com += points[i]
  com /= n
  return com

def add_vector_to_points(vector, points):
  new_points = []
  n = len(points)
  for i in range(n):
    new_points.append(vector+points[i])
  return new_points

def subtract_vector_from_points(vector, points):
  new_points = []
  n = len(points)
  for i in range(n):
    new_points.append(points[i]-vector)
  return new_points

def find_rmsd_anchor_from_vectors(points1, points2, vector, radius, TOL=0.0001):
  '''Given a set of points in the receptor (points1) and ligand (points2)
  the algorithm uses a Babylonian algorithm to choose a location for points2
  that lies along the vector
  by summing them until they cross the radius. Then, the ligand will be placed
  at that location - the anchor.
  Input:
   - origin: 3-membered numpy array representing the center of the spherical 
  milestone
   - radius: float representing the radius of the milestone
   - vectors: a list of 3-membered arrays representing vectors that define the
   'pathway' out of the receptor's active site
  Output:
   - anchor: 3-array representing the location of the anchor.
  '''
  max_counter = 10000
  counter = 0
  assert len(points1) == len(points2), "there must be the same number of points in the points1 and points2 variables."
  assert np.linalg.norm(vector) > 0.0, "vector must not be length zero"
  n = len(points1)
  com1 = get_points_com(points1)
  com2 = get_points_com(points2)
  points_diff_xtal = get_points_diff(points1, points2)
  xtal_rmsd = points_rmsd(points_diff_xtal)
  points1_origin = subtract_vector_from_points(com1, points1)
  points2_origin = subtract_vector_from_points(com2, points2)
  points_diff_origin = get_points_diff(points1_origin, points2_origin)
  origin_rmsd = points_rmsd(points_diff_origin)
  assert radius > origin_rmsd, "Specified RMSD milestone radius impossibly small: no way to obtain specified radius. Radius: %f, RMSD: %f" % (radius, origin_rmsd)
  points2_vector = add_vector_to_points(vector, points2_origin)
  vec_lower = np.array([0.0, 0.0, 0.0])
  vec_higher = deepcopy(vector)
  
  while True:
    points2_vec_higher = add_vector_to_points(vec_higher, points2_origin)
    points_diff = get_points_diff(points1_origin, points2_vec_higher)
    rmsd = points_rmsd(points_diff)
    #print "rmsd:", rmsd
    if (rmsd - radius)**2 < TOL**2: # then we have found a radius that works
      #new_points2 = add_vector_to_points(com1, add_vector_to_points(vec_higher, points2_origin))
      anchor = com1 + vec_higher
      break
    elif radius > rmsd:
      vector *= 2
      vec_higher = deepcopy(vector)
    else:
      vector = 0.5*(vec_higher + vec_lower)
      points2_vector = add_vector_to_points(vector, points2_origin)
      points_diff_vector = get_points_diff(points1_origin, points2_vector)
      rmsd2 = points_rmsd(points_diff_vector)
      if radius > rmsd2: # need to go higher
        vec_lower = deepcopy(vector)
      else: # need to go lower
        vec_higher = deepcopy(vector)
  
    counter += 1
    assert counter < max_counter, "unable to converge on proper radius in RMSD milestone. Consider using larger TOL value."
  
  return com1, com2, anchor

def find_planar_z_anchor_from_vectors(origin, offset, vectors):
  '''Given an origin and radius of the milestone, will follow the vectors list
  by summing them until they cross the radius. Then, the ligand will be placed
  at that location - the anchor.
  Input:
   - origin: 3-membered numpy array representing the center of the spherical 
  milestone
   - offset: float representing the offset of the milestone
   - vectors: a list of 3-membered arrays representing vectors that define the
   'pathway' out of the receptor's active site
  Output:
   - anchor: 3-array representing the location of the anchor.
  '''
  summed_vectors = np.array([0.0,0.0,0.0]) # starting from the very beginning of the origin
  for i in range(len(vectors)):
    cur_vector = vectors[i]
    #cur_vector_prime = cur_vector / np.linalg.norm(cur_vector) # TODO: remove
    if i == len(vectors)-1 or summed_vectors[2] + cur_vector[2] > offset: # then this is the last vector, so the anchor has to be placed in line with this one
      '''# TODO: remove
      a = np.dot(cur_vector_prime,cur_vector_prime) # 'a' for the quadratic formula
      b = 2.0 * np.dot(summed_vectors,cur_vector_prime) # 'b' for the quadratic formula
      c = np.dot(summed_vectors,summed_vectors) - radius**2 # 'c' for the quadratic formula
      roots = quadratic_formula(a,b,c)
      '''
      a = offset - summed_vectors[2] - origin[2]
      assert cur_vector[2] != 0.0, "Must not have a placement vector parallel to the milestones."
      cur_vector_znorm = cur_vector / cur_vector[2]
      anchor = origin + summed_vectors + cur_vector_znorm * a # get the root with the smallest absolute value from the quadratic formula
      return anchor
    summed_vectors += cur_vector
  return anchor
    
def find_ellipsoidal_anchor_from_vectors(locus_1, locus_2, nu):
  '''Given the locations of the two loci of the milestone, the ligand will be 
  placed at the ellipsoidal coordinate 'nu' between the two points.
  Input:
   - locus_1: 3-membered numpy array representing the center of locus 1 
   - locus_2: 3-membered numpy array representing the center of locus 2
   - nu: the nu coordinate
  Output:
   - anchor: 3-array representing the location of the anchor.
  '''
  diff_vec = locus_2 - locus_1
  anchor = diff_vec * nu
  return anchor

'''
def find_RMSD_anchor_from_vectors(origin, radius, vectors):
  ''Given an origin and radius of the milestone, will follow the vectors list
  by summing them until they cross the radius. Then, the ligand will be placed
  at that location - the anchor.
  Input:
   - origin: 3-membered numpy array representing the center of the spherical 
  milestone
   - radius: float representing the radius of the milestone
   - vectors: a list of 3-membered arrays representing vectors that define the
   'pathway' out of the receptor's active site
  Output:
   - anchor: 3-array representing the location of the anchor.
  ''
  summed_vectors = np.array([0.0,0.0,0.0]) # starting from the very beginning of the origin 
  for i in range(len(vectors)):
    cur_vector = vectors[i]
    cur_vector_prime = cur_vector / np.linalg.norm(cur_vector)
    if i == len(vectors)-1 or np.linalg.norm(summed_vectors + cur_vector) > radius: # then this is the last vector, so the anchor has to be placed in line with this one
      a = np.dot(cur_vector_prime,cur_vector_prime) # 'a' for the quadratic formula
      b = 2.0 * np.dot(summed_vectors,cur_vector_prime) # 'b' for the quadratic formula
      c = np.dot(summed_vectors,summed_vectors) - radius**2 # 'c' for the quadratic formula
      roots = quadratic_formula(a,b,c)
      anchor = origin + summed_vectors + cur_vector_prime * min(map(abs,roots)) # get the root with the smallest absolute value from the quadratic formula
      return anchor
    summed_vectors += cur_vector
  return anchor
''' # remove
    
def generate_spherical_milestones(seekrcalc, atom_indices, origin, radius_list, siteid, vectors, absolute=False, lig_center=None):
  '''Create a list of concentric spherical milestones around an atom selection.
  Input:
   - atom_indices: a list of integers representing the atom indices that this milestone is centered on
   - origin: 3-array that should be equal (or approximately equal) to the coordinates of the center of mass of the atom indices in the input structure
   - radius_list: a list of floats representing the radii of the concentric spheres in the list of milestones
   - siteid: integer representing the index of this list of milestones
   - vectors: a list of lists of floats, or of numpy arrays, that define the pathway out of the active site, and in turn, the locations where the ligand will be placed
   - absolute: boolean defining whether this milestone will remain stationary in space regardless of how the receptor moves  
  Output:
   - milestones: a list of milestone objects
  '''
    
  milestones = []
  lowest_radius = radius_list[0] # the lowest radius
  spheres_in_site = len(radius_list)
  for i in range(spheres_in_site): # loop through all the milestone radii
    radius = radius_list[i]
    anchor = find_spherical_anchor_from_vectors(origin, radius, vectors)# now we need to find the anchor - the location where the ligand will be placed
    anchor_tuple = tuple(anchor) # TODO: remove?
    milestone = Concentric_Spherical_Milestone(i, siteid, absolute, md=seekrcalc.project.md, bd=seekrcalc.project.bd)
    milestone.center_atom_indices = atom_indices
    milestone.center_vec = origin
    milestone.lig_center = lig_center
    milestone.radius = radius
    milestone.anchor = anchor
    if i > 0:
      milestone.neighbors.append(i-1) # if we are the furthest milestone in either direction, then include no neighbors
    #else:
      #if not k_off: milestone.end = True
    if i < spheres_in_site - 1:
      milestone.neighbors.append(i+1)
      milestone.bd = False
    else: # then this is the outermost one, so set BD to true
      milestone.bd = True
      milestone.bd_adjacent = milestones[-1] # make the adjacent milestone the previously created md milestone
      milestone.end = True
      milestone.md = False
    milestone.fullname = "%d_%d_%s_%.1f_%.1f_%.1f" % (i+(spheres_in_site*int(siteid)), i, siteid, anchor[0], anchor[1], anchor[2])
    milestones.append(milestone)
  return milestones

def print_spherical_milestone_info(milestones):
  '''Print the information in all milestones for debug and verbosity purposes.
  Input:
   - milestones: a list of milestone objects to print information about.
  Output:
   - None
   '''
  for milestone in milestones:
    print "Milestone index:", milestone.index, "siteid:", milestone.siteid, "origin:", milestone.center_vec, "center_atoms:", milestone.center_atom_indices
    print "  num neighbors:", len(milestone.neighbors), "anchor:", milestone.anchor, "md:", milestone.md, "bd:", milestone.bd, "bd_adjacent:", milestone.bd_adjacent
    print "  end:", milestone.end, "fullname:", milestone.fullname, "radius:", milestone.radius
  return

def generate_planar_z_milestones(seekrcalc, atom_indices, origin, offset_list, 
                                 siteid, vectors=[np.array([0.0, 0.0, 1.0])], 
                                 absolute=False, lig_center=None):

  '''Create a list of stacked planar milestones around an atom selection.
  Input:
   - atom_indices: a list of integers representing the atom indices that this 
     milestone is centered on
   - origin: 3-array that should be equal (or approximately equal) to the 
     coordinates of the center of mass of the atom indices in the input 
     structure
   - offset_list: a list of floats representing the offsets of the planes in 
     the list of milestones
   - siteid: integer representing the index of this list of milestones
   - vectors: a list of lists of floats, or of numpy arrays, that define the 
     pathway out of the active site, and in turn, the locations where the 
     ligand will be placed
   - absolute: boolean defining whether this milestone will remain stationary 
     in space regardless of how the receptor moves  
  Output:
   - milestones: a list of milestone objects
  '''
    
  milestones = []
  #lowest_radius = radius_list[0] # the lowest radius #TODO: remove
  #total_spherical_milestones = 0 #TODO: remove
  num_planes = len(offset_list)
  #cur_vector = vectors[0] # TODO: remove
  for i in range(num_planes): # loop through all the planar milestones
    offset = offset_list[i]
    #anchor = find_spherical_anchor_from_vectors(origin, radius, vectors)# now we need to find the anchor - the location where the ligand will be placed # TODO: marked for removal
    anchor = find_planar_z_anchor_from_vectors(origin, offset, vectors)
    #anchor_tuple = tuple(anchor) # TODO: remove?
    milestone = Planar_Z_Milestone(i, siteid, absolute, md=seekrcalc.project.md, bd=seekrcalc.project.bd)
    milestone.center_atom_indices = atom_indices
    milestone.center_vec = origin
    milestone.lig_center = lig_center
    milestone.offset = offset
    milestone.anchor = anchor
    milestone.bd = False
    milestone.md = True
    
    if i > 0:
      milestone.neighbors.append(i-1) # if we are the furthest milestone in either direction, then include no neighbors
    
    if i < num_planes - 1:
      milestone.neighbors.append(i+1)
    
    milestone.fullname = "%d_%d_%s_%.1f_%.1f_%.1f" % (i+(num_planes*int(siteid)), i, siteid, anchor[0], anchor[1], anchor[2])
    milestones.append(milestone)
  return milestones

def print_planar_z_milestone_info(milestones):
  '''Print the information in all milestones for debug and verbosity purposes.
  Input:
   - milestones: a list of milestone objects to print information about.
  Output:
   - None
   '''
  for milestone in milestones:
    print "Milestone index:", milestone.index, "siteid:", milestone.siteid, \
        "origin:", milestone.center_vec, "center_atoms:", \
        milestone.center_atom_indices
    print "  num neighbors:", len(milestone.neighbors), "anchor:", \
        milestone.anchor, "md:", milestone.md, "bd:", milestone.bd, \
        "bd_adjacent:", milestone.bd_adjacent
    print "  end:", milestone.end, "fullname:", milestone.fullname, "offset:", \
        milestone.offset
  return

def generate_ellipsoidal_milestones(seekrcalc, atom_indices_locus_1, 
                                    atom_indices_locus_2, origin_1, origin_2, 
                                    nu_list, siteid, absolute=False, 
                                    lig_center=None):
  '''Create a list of concentric spherical milestones around an atom selection.
  Input: UPDATE
   - atom_indices: a list of integers representing the atom indices that this milestone is centered on
   - origin: 3-array that should be equal (or approximately equal) to the coordinates of the center of mass of the atom indices in the input structure
   - radius_list: a list of floats representing the radii of the concentric spheres in the list of milestones
   - siteid: integer representing the index of this list of milestones
   - vectors: a list of lists of floats, or of numpy arrays, that define the pathway out of the active site, and in turn, the locations where the ligand will be placed
   - absolute: boolean defining whether this milestone will remain stationary in space regardless of how the receptor moves  
  Output:
   - milestones: a list of milestone objects
  '''
    
  milestones = []
  nus_in_site = len(nu_list)
  for i in range(nus_in_site): # loop through all the milestone nu's
    nu = nu_list[i]
    anchor = find_ellipsoidal_anchor_from_vectors(origin_1, origin_2, nu)# now we need to find the anchor - the location where the ligand will be placed
    milestone = Ellipsoidal_Milestone(i, siteid, absolute, md=seekrcalc.project.md, bd=seekrcalc.project.bd)
    milestone.locus_1_atom_indices = atom_indices_locus_1
    milestone.locus_1_vec = origin_1
    milestone.locus_2_atom_indices = atom_indices_locus_2
    milestone.locus_2_vec = origin_2
    milestone.lig_center = lig_center
    milestone.nu = nu
    milestone.anchor = anchor
    milestone.bd = False
    if i > 0:
      milestone.neighbors.append(i-1) # if we are the furthest milestone in either direction, then include no neighbors
    
    if i < nus_in_site - 1:
      milestone.neighbors.append(i+1)
      
    milestone.fullname = "%d_%d_%s_%.1f_%.1f_%.1f" % (i+(nus_in_site*int(siteid)), i, siteid, anchor[0], anchor[1], anchor[2])
    milestones.append(milestone)
  return milestones

def print_ellipsoidal_milestone_info(milestones):
  '''Print the information in all milestones for debug and verbosity purposes.
  Input:
   - milestones: a list of milestone objects to print information about.
  Output:
   - None
   '''
  for milestone in milestones:
    print "Milestone index:", milestone.index, "siteid:", milestone.siteid
    print "origin_1:", milestone.locus_1_vec, "center_atoms_locus_1:", milestone.locus_1_atom_indices
    print "origin_2:", milestone.locus_2_vec, "center_atoms_locus_2:", milestone.locus_2_atom_indices
    print "num neighbors:", len(milestone.neighbors), "anchor:", milestone.anchor
    print "md:", milestone.md, "bd:", milestone.bd, "bd_adjacent:", milestone.bd_adjacent
    print "end:", milestone.end, "fullname:", milestone.fullname, "nu:", milestone.nu
  return

def generate_RMSD_milestones(seekrcalc, atom_indices1, atom_indices2, 
                             atom_points1, atom_points2, radius_list, 
                             siteid, vector, absolute=False, 
                             lig_center=None):
  '''Create a list of concentric RMSD milestones around an atom selection.
  Input:
   - atom_indices: a list of integers representing the atom indices that this milestone is centered on
   - origin: 3-array that should be equal (or approximately equal) to the coordinates of the center of mass of the atom indices in the input structure
   - radius_list: a list of floats representing the radii of the concentric RMSD hyperspheres in the list of milestones
   - siteid: integer representing the index of this list of milestones
   - vectors: a list of lists of floats, or of numpy arrays, that define the pathway out of the active site, and in turn, the locations where the ligand will be placed
   - absolute: boolean defining whether this milestone will remain stationary in space regardless of how the receptor moves  
  Output:
   - milestones: a list of milestone objects
  '''
    
  milestones = []
  lowest_radius = radius_list[0] # the lowest radius
  #total_spherical_milestones = 0 #TODO: remove
  hyperspheres_in_site = len(radius_list)
  #cur_vector = vectors[0] # TODO: remove
  for i in range(hyperspheres_in_site): # loop through all the milestone radii
    radius = radius_list[i]
    com1, com2, anchor = find_rmsd_anchor_from_vectors(atom_points1, atom_points2, vector, radius)# now we need to find the anchor - the location where the ligand will be placed
    milestone = RMSD_Milestone(i, siteid, absolute, md=seekrcalc.project.md, bd=seekrcalc.project.bd)
    milestone.atom_indices1 = atom_indices1
    milestone.atom_points1 = atom_points1
    milestone.atom_vec1 = com1
    milestone.atom_indices2 = deepcopy(atom_indices2)
    milestone.atom_points2 = atom_points2
    milestone.atom_vec2 = anchor
    milestone.radius = radius
    milestone.anchor = anchor
    milestone.lig_center = com2
    if i > 0:
      milestone.neighbors.append(i-1) # if we are the furthest milestone in either direction, then include no neighbors
    #else:
      #if not k_off: milestone.end = True
    if i < hyperspheres_in_site - 1:
      milestone.neighbors.append(i+1)
      milestone.bd = False
    else: # then this is the outermost one, so set BD to true
      milestone.bd = True
      milestone.bd_adjacent = milestones[-1] # make the adjacent milestone the previously created md milestone
      milestone.end = True
      milestone.md = False
    milestone.fullname = "%d_%d_%s_%.1f_%.1f_%.1f" % (i+(hyperspheres_in_site*int(siteid)), i, siteid, anchor[0], anchor[1], anchor[2])
    milestones.append(milestone)
    
  return milestones

def print_RMSD_milestone_info(milestones):
  '''Print the information in all milestones for debug and verbosity purposes.
  Input:
   - milestones: a list of milestone objects to print information about.
  Output:
   - None
   '''
  for milestone in milestones:
    print "Milestone index:", milestone.index, "siteid:", milestone.siteid
    print "atom_vec1:", milestone.atom_vec1, "atom_indices1:", milestone.atom_indices1
    print "atom_vec2:", milestone.atom_vec2, "atom_indices2:", milestone.atom_indices2
    print "  num neighbors:", len(milestone.neighbors), "anchor:", milestone.anchor
    print "md:", milestone.md, "bd:", milestone.bd, "bd_adjacent:", milestone.bd_adjacent
    print "  end:", milestone.end, "fullname:", milestone.fullname, "radius:", milestone.radius
  return


class Test_milestones(unittest.TestCase):
  # several test cases to ensure the functions in this module are working properly
  def test_main(self): # test this function
    print "WARNING: this module does not have comprehensive unit tests"
    
  def test_quadratic_formula(self):
    self.assertEqual(quadratic_formula(1.0,4.0,3.0), [-3.0,-1.0]) # test negative roots
    self.assertEqual(quadratic_formula(3.0,-13.0,12.0), [4.0/3.0,3.0]) # test positive roots
    self.assertEqual(quadratic_formula(1.0,-2.0,0.0), [0.0,2.0]) # test zero roots
    self.assertEqual(quadratic_formula(1.0,0.0,0.0), [0.0,0.0]) # test zero roots
    self.assertEqual(quadratic_formula(1.0,0.0,4.0), [complex(0,-2.0),complex(0,2.0)]) # test imaginary roots
    return

  def test_find_rad_of_skewed_vectors(self):
    A = np.array([0.0,1.0,0.0])
    B = np.array([1.0,0.0,0.0])
    B2 = np.array([1.0,1.0,0.0])
    origin = np.array([0,0,0])
    ''' # TODO: complete this later
    self.assertEqual(find_spherical_anchor_from_vectors(origin, A, 2.0), 1.0) # check a vector against itself
    self.assertEqual(find_spherical_anchor_from_vectors(A, B, 2.0), sqrt(3.0)) # check a vector against right angle
    self.assertEqual(find_spherical_anchor_from_vectors(2.0*A, B*2.5, 3.0), sqrt(5.0)) # check a vector against right angle
    self.assertEqual(find_spherical_anchor_from_vectors(B2, A, 2.0), sqrt(3.0)-1.0) # check a vector against 45 degree angle
    '''
    return
  
  def test_get_points_diff(self):
    '''Find vectors (points_diff) that points from points1 to points2'''
    #print "get_points_diff(np.array([5.0, 5.0, 5.0]), np.array([2.0, 3.0, 4.0])):", get_points_diff(np.array([5.0, 5.0, 5.0]), np.array([2.0, 3.0, 4.0]))
    self.assertTrue(np.array_equal([np.array([3.0, 2.0, 1.0])], get_points_diff([np.array([5.0, 5.0, 5.0])], [np.array([2.0, 3.0, 4.0])])))
  
  def test_points_rmsd(self):
    '''Finds the RMSD between two sets of points'''
    #print "points_rmsd(np.array([3.0, 4.0, 5.0])):", points_rmsd(np.array([3.0, 4.0, 5.0]))
    self.assertAlmostEqual(np.sqrt(50.0), points_rmsd(np.array([3.0, 4.0, 5.0])))
    

  def test_get_points_com(self):
    '''compute the center of mass of a vector.'''
  
    points = [np.array([-1.0, 2.0, 0.0]), np.array([1.0, 2.0, 0.0]),
            np.array([0.0, 1.0, 0.0]), np.array([0.0, 3.0, 0.0]),
            np.array([0.0, 2.0, -1.0]), np.array([0.0, 2.0, 1.0])]
    #print "get_points_com(points):", get_points_com(points)
    self.assertTrue(np.array_equal(np.array([0.0, 2.0, 0.0]), get_points_com(points)))
  

  def test_add_vector_to_points(self):
    vector = np.array([-3.0, 4.0, 5.0])
    points = [np.array([-1.0, 1.0, 0.0]), np.array([1.0, 1.0, 0.0])]
    result = [np.array([-4.0, 5.0, 5.0]), np.array([-2.0, 5.0, 5.0])]
    #print "result:", result
    #print "add_vector_to_points(vector, points):", add_vector_to_points(vector, points)
    self.assertTrue(np.array_equal(result[0], add_vector_to_points(vector, points)[0]))
    self.assertTrue(np.array_equal(result[1], add_vector_to_points(vector, points)[1]))

  def subtract_vector_from_points(self):
    vector = np.array([-3.0, 4.0, 5.0])
    points = [np.array([-1.0, 1.0, 0.0]), np.array([1.0, 1.0, 0.0])]
    result = [np.array([2.0, -3.0, -5.0]), np.array([4.0, -3.0, -5.0])]
    #print "result:", result
    #print "add_vector_to_points(vector, points):", add_vector_to_points(vector, points)
    self.assertEqual(result, add_vector_to_points(vector, points))

if __name__ == "__main__":
  unittest.main()
