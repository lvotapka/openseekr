'''
Created on May 9, 2018

@author: lvotapka

Procedures for creating milestone objects
'''

import numpy as np
from math import sqrt
import cmath
import unittest
import xml.etree.ElementTree as ET
from base import strBool
from simtk.unit import nanometer, Quantity



verbose = True

class Milestone_System():
    '''An object for system information stored for each milestone.'''
    def __init__(self):
        self.wet_holo_pdb_filename = ''
        self.dry_holo_pdb_filename = ''
        self.prmtop_filename = ''
        self.inpcrd_filename = ''
        self.umbrella_pdb_filename = ''
        self.system = None
        self.simulation = None
        return
    
    def serialize(self, xmlMilestone_system):
        xmlWet_holo_pdb_filename = ET.SubElement(xmlMilestone_system, 'wet_holo_pdb_filename')
        xmlWet_holo_pdb_filename.text = str(self.wet_holo_pdb_filename)
        xmlDry_holo_pdb_filename = ET.SubElement(xmlMilestone_system, 'dry_holo_pdb_filename')
        xmlDry_holo_pdb_filename.text = str(self.dry_holo_pdb_filename)
        xmlPrmtop_filename = ET.SubElement(xmlMilestone_system, 'prmtop_filename')
        xmlPrmtop_filename.text = str(self.prmtop_filename)
        xmlInpcrd_filename = ET.SubElement(xmlMilestone_system, 'inpcrd_filename')
        xmlInpcrd_filename.text = str(self.inpcrd_filename)
        xmlUmbrella_pdb_filename = ET.SubElement(xmlMilestone_system, 'umbrella_pdb_filename')
        xmlUmbrella_pdb_filename.text = str(self.umbrella_pdb_filename)
        return
    
    def deserialize(self, xmlTree):
        self.wet_holo_pdb_filename = xmlTree.find('wet_holo_pdb_filename').text
        self.dry_holo_pdb_filename = xmlTree.find('dry_holo_pdb_filename').text
        self.prmtop_filename = xmlTree.find('prmtop_filename').text
        self.inpcrd_filename = xmlTree.find('inpcrd_filename').text
        self.umbrella_pdb_filename = xmlTree.find('umbrella_pdb_filename').text
        return

class Milestone():
    '''Milestone superclass. It represents a surface in phase space that is monitored for crossings'''
    def __init__(self, anchor, dimensions, index, siteid):
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
    def __init__(self, index=0, siteid=0, absolute='False', md=True, bd=False):
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
        self.radius = 0.0
        self.wet_holo_filename = ''
        self.dry_holo_filename = ''
        self.atom_selection_1 = None
        self.atom_selection_2 = None
        self.openmm = Milestone_System()
        self.config = None
        self.box_vectors = None
        return
        
    def serialize(self, xmlMilestone):
        xmlType = ET.SubElement(xmlMilestone, 'type')
        xmlType.text = 'concentric_spherical'
        xmlFullname = ET.SubElement(xmlMilestone, 'fullname')
        xmlFullname.text = str(self.fullname)
        xmlDirectory = ET.SubElement(xmlMilestone, 'directory')
        xmlDirectory.text = str(self.directory)
        xmlAnchor = ET.SubElement(xmlMilestone, 'anchor')
        xmlAnchor.text = str(self.anchor)
        xmlNeighbors = ET.SubElement(xmlMilestone, 'neighbors')
        for neighbor in self.neighbors:
            xmlNeighbor = ET.SubElement(xmlNeighbors, 'neighbor')
            xmlNeighbor.text = str(neighbor)
        xmlIndex = ET.SubElement(xmlMilestone, 'index')
        xmlIndex.text = str(self.index)
        xmlSiteID = ET.SubElement(xmlMilestone, 'siteid')
        xmlSiteID.text = str(self.siteid)
        xmlAbsolute = ET.SubElement(xmlMilestone, 'absolute')
        xmlAbsolute.text = str(self.absolute)
        xmlMd = ET.SubElement(xmlMilestone, 'md')
        xmlMd.text = str(self.md)
        xmlBd = ET.SubElement(xmlMilestone, 'bd')
        xmlBd.text = str(self.bd)
        xmlBd_adjacent = ET.SubElement(xmlMilestone, 'bd_adjacent')
        if self.bd_adjacent:
            xmlBd_adjacent.text = str(self.bd_adjacent.index)
        else:
            xmlBd_adjacent.text = ''
        
        xmlEnd = ET.SubElement(xmlMilestone, 'end')
        xmlEnd.text = str(self.end)
        xmlCenter_atom_indices = ET.SubElement(xmlMilestone, 'center_atom_indices')
        xmlCenter_atom_indices.text = ', '.join(list(map(str, self.center_atom_indices)))
        xmlCenter_vec = ET.SubElement(xmlMilestone, 'center_vec')
        xmlCenter_vec.text = ', '.join(list(map(str, self.center_vec)))
        xmlRadius = ET.SubElement(xmlMilestone, 'radius')
        xmlRadius.text = str(self.radius)
        xmlWet_holo_pdb_filename = ET.SubElement(xmlMilestone, 'wet_holo_filename')
        xmlWet_holo_pdb_filename.text = str(self.wet_holo_filename)
        xmlDry_holo_pdb_filename = ET.SubElement(xmlMilestone, 'dry_holo_filename')
        xmlDry_holo_pdb_filename.text = str(self.dry_holo_filename)
        xmlAtom_selection_1 = ET.SubElement(xmlMilestone, 'atom_selection_1')
        xmlAtom_selection_1.text = ', '.join(list(map(str, self.atom_selection_1)))
        xmlAtom_selection_2 = ET.SubElement(xmlMilestone, 'atom_selection_2')
        xmlAtom_selection_2.text = ', '.join(list(map(str, self.atom_selection_2)))
        xmlOpenMM = ET.SubElement(xmlMilestone, 'openmm')
        xmlOpenMM.text = self.openmm.serialize(xmlOpenMM)
        #xmlConfig = ET.SubElement(xmlMilestone, 'config')
        #xmlConfig.text = str(self.config)
        xmlBox_vectors = ET.SubElement(xmlMilestone, 'box_vectors')
        if self.box_vectors is not None:
            box_vectors_unitless = self.box_vectors.value_in_unit(nanometer)
            xmlA = ET.SubElement(xmlBox_vectors, 'A')
            xmlAx = ET.SubElement(xmlA, 'x')
            xmlAx.text = str(box_vectors_unitless[0][0])
            xmlAy = ET.SubElement(xmlA, 'y')
            xmlAy.text = str(box_vectors_unitless[0][1])
            xmlAz = ET.SubElement(xmlA, 'z')
            xmlAz.text = str(box_vectors_unitless[0][2])
            xmlB = ET.SubElement(xmlBox_vectors, 'B')
            xmlBx = ET.SubElement(xmlB, 'x')
            xmlBx.text = str(box_vectors_unitless[1][0])
            xmlBy = ET.SubElement(xmlB, 'y')
            xmlBy.text = str(box_vectors_unitless[1][1])
            xmlBz = ET.SubElement(xmlB, 'z')
            xmlBz.text = str(box_vectors_unitless[1][2])
            xmlC = ET.SubElement(xmlBox_vectors, 'C')
            xmlCx = ET.SubElement(xmlC, 'x')
            xmlCx.text = str(box_vectors_unitless[2][0])
            xmlCy = ET.SubElement(xmlC, 'y')
            xmlCy.text = str(box_vectors_unitless[2][1])
            xmlCz = ET.SubElement(xmlC, 'z')
            xmlCz.text = str(box_vectors_unitless[2][2])
        else:
            xmlBox_vectors.text = ''
        return
    
    def deserialize(self, xmlTree):
        self.fullname = xmlTree.find('fullname').text
        self.directory = xmlTree.find('directory').text
        self.anchor = xmlTree.find('anchor').text
        self.index = int(xmlTree.find('index').text)
        self.siteid = int(xmlTree.find('siteid').text)
        self.absolute = xmlTree.find('absolute').text
        self.md = strBool(xmlTree.find('md').text)
        self.bd = strBool(xmlTree.find('bd').text)
        # temporarily store index, fill with object later
        bd_adjacent_str = xmlTree.find('bd_adjacent').text
        if bd_adjacent_str:
            self.bd_adjacent = int(bd_adjacent_str) 
        self.end = strBool(xmlTree.find('end').text)
        center_atom_indices_str = xmlTree.find('center_atom_indices').text
        if center_atom_indices_str:
            center_atom_indices_str = center_atom_indices_str.replace(' ', '')
            self.center_atom_indices = list(
                map(int, center_atom_indices_str.split(',')))
        center_vec_str = xmlTree.find('center_vec').text
        if center_vec_str:
            center_vec_str = center_vec_str.replace(' ', '')
            self.center_vec = list(
                map(float, center_vec_str.split(',')))
        self.radius = float(xmlTree.find('radius').text)
        self.wet_holo_filename = xmlTree.find('wet_holo_filename').text
        self.dry_holo_filename = xmlTree.find('dry_holo_filename').text
        atom_selection_1_str = xmlTree.find('atom_selection_1').text
        if atom_selection_1_str:
            atom_selection_1_str = atom_selection_1_str.replace(' ', '')
            self.atom_selection_1 = list(
                map(int, atom_selection_1_str.split(',')))
        atom_selection_2_str = xmlTree.find('atom_selection_2').text
        if atom_selection_2_str:
            atom_selection_2_str = atom_selection_2_str.replace(' ', '')
            self.atom_selection_2 = list(
                map(int, atom_selection_2_str.split(',')))
        self.openmm.deserialize(xmlTree.find('openmm'))
        xmlBox_vectors = xmlTree.find('box_vectors')
        if xmlBox_vectors.text is not None:
            xmlA = xmlBox_vectors.find('A')
            xmlAx = float(xmlA.find('x').text)
            xmlAy = float(xmlA.find('y').text)
            xmlAz = float(xmlA.find('z').text)
            xmlB = xmlBox_vectors.find('B')
            xmlBx = float(xmlB.find('x').text)
            xmlBy = float(xmlB.find('y').text)
            xmlBz = float(xmlB.find('z').text)
            xmlC = xmlBox_vectors.find('C')
            xmlCx = float(xmlC.find('x').text)
            xmlCy = float(xmlC.find('y').text)
            xmlCz = float(xmlC.find('z').text)
            self.box_vectors = Quantity([[xmlAx, xmlAy, xmlAz], 
                                         [xmlBx, xmlBy, xmlBz],
                                         [xmlCx, xmlCy, xmlCz]], 
                                         unit=nanometer)
        else:
            self.box_vectors = None
        return

def deserialize_milestones(xmlTree):
    '''
    Deserialize a milestone section of the SeekrCalculation XML file
    '''
    milestones_list = []
    for milestone_xml in xmlTree:
        milestone_type = milestone_xml.find('type').text
        if milestone_type == None or \
                milestone_type == 'concentric_spherical':
            milestone = Concentric_Spherical_Milestone()
            milestone.deserialize(milestone_xml)
        else:
            raise Exception("milestone type not implemented: %s" % \
                            milestone_type)
            
        milestones_list.append(milestone)
            
    for milestone in milestones_list:
        bd_adjacent_index = milestone.bd_adjacent
        if bd_adjacent_index is not None:
            milestone.bd_adjacent = milestones_list[bd_adjacent_index]
            
    return milestones_list

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

def find_anchor_from_vectors(origin, radius, vectors):
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
    summed_vectors = np.array([0.0,0.0,0.0]) # starting from the very beginning of the
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

def generate_spherical_milestones(seekrcalc, atom_indices, origin, radius_list, siteid, vectors, absolute=False):
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
    #total_spherical_milestones = 0 #TODO: remove
    spheres_in_site = len(radius_list)
    #cur_vector = vectors[0] # TODO: remove
    for i in range(spheres_in_site): # loop through all the milestone radii
        radius = radius_list[i]
        anchor = find_anchor_from_vectors(origin, radius, vectors)# now we need to find the anchor - the location where the ligand will be placed
        anchor_tuple = tuple(anchor)
        milestone = Concentric_Spherical_Milestone(i, siteid, absolute, md=seekrcalc.project.md, bd=seekrcalc.project.bd)
        milestone.center_atom_indices = atom_indices
        milestone.center_vec = origin
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
        print "Milestone index:", milestone.index, "siteid:", milestone.siteid, "origin:", milestone.center_vec, "center_atoms:", milestone.center_atom_indices,
        print "  num neighbors:", len(milestone.neighbors), "anchor:", milestone.anchor, "md:", milestone.md, "bd:", milestone.bd, "bd_adjacent:", milestone.bd_adjacent
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
        self.assertEqual(find_anchor_from_vectors(origin, A, 2.0), 1.0) # check a vector against itself
        self.assertEqual(find_anchor_from_vectors(A, B, 2.0), sqrt(3.0)) # check a vector against right angle
        self.assertEqual(find_anchor_from_vectors(2.0*A, B*2.5, 3.0), sqrt(5.0)) # check a vector against right angle
        self.assertEqual(find_anchor_from_vectors(B2, A, 2.0), sqrt(3.0)-1.0) # check a vector against 45 degree angle
        '''
        return
