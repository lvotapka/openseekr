#!/usr/bin/python

"""
bd.py

creates the necessary files to run a BD simulation using Browndye


"""
import sys
import os 
import math
import shutil 
import subprocess 
import glob #, make_fxd
import numpy as np

import pdb2 as pdb
# needed to keep track of separate structure objects
from copy import deepcopy
import unittest
import re
import apbs
from xml.dom.minidom import Document
import copy
import random
import xml.etree.cElementTree as ET # for writing xml files
from xml.dom import minidom
from adv_template import Adv_template, File_template


verbose = True

parser = pdb.Big_PDBParser()


#browndye_bin = '/soft/browndye/latest/bin/'
RXN_FILENAME = 'rxns.xml'
INPUT_FILENAME = 'input.xml'
DEFAULT_TEMP = 298.0
empty_pqrxml = "./empty.pqrxml"
#FHPD_NUMTRAJ = 1000 # should be set within the settings file

"""
 inputs needed:
 -PQR files

 steps needed to run a BD simulation with Browndye
 1. Run Electrostatics -> get E. potential & debye length
 2. Determine reaction criteria

"""

empty_pqrxml_template = """<roottag>
</roottag>"""

extract_bd_frames_template = """#!/usr/bin/python

\"\"\"
Given a series of trajectories, will extract out a pqr file for each 
encounter event to get a First Hitting Point Distribution.
\"\"\"

import os, sys, re

# CONSTANTS -- anything in this section can be 
# safely modified by automation or manually

# the directory that contains all the FHPD BD trajectories
trajdir = "$TRAJDIR"
workdir = "$WORKDIR"
pqrxml0 = os.path.join(trajdir,"$PQRXML0") # pqrxml file of the receptor
pqrxml1 = os.path.join(trajdir,"$PQRXML1") # pqrxml file of the ligand
empty = "$EMPTY"
# the name of the site in the b_surface BD rxns file
sitename = "$SITENAME"
number_of_trajs = $NUMBER_OF_TRAJS
max_structures = 1000
counter = 0

# SCRIPT -- Nothing beyond this point should 
# need to be modified (unless there's a bug)
def modify_pqr(pqr_filename):
  '''edits the pqr file by removing any line except lines that begin with ATOM because BrownDye has weird pqr output'''
  pqrread = open(pqr_filename, 'r')
  raw_pqr = pqrread.readlines()
  pqrread.close
  pqrwrite = open(pqr_filename, 'w')
  for line in raw_pqr:
    if re.match("ATOM",line):
      pqrwrite.write(line)
  pqrwrite.close()
  return

print "now extracting all successful reaction numbers..."
# create the working directory
if not os.path.exists(workdir): os.mkdir(workdir)
quitting = False
for i in range(number_of_trajs):
  if quitting: break
  print "extracting trajectories from traj number:", i
  outputfilename = os.path.join(workdir, "rxn_output%d.txt" % i)
  traj_filename = os.path.join(trajdir, "traj%d.xml" % i)
  trajindex_filename = os.path.join(trajdir, "traj%d.index.xml" % i)
  cmd = "echo 'Browndye Trajectory number' > %s; \
      process_trajectories -traj %s -index %s -srxn %s >> %s" \
      % (outputfilename, traj_filename, 
         trajindex_filename, sitename, outputfilename)
  os.system(cmd) # run the command to extract successful trajectories
  print "running command:", cmd
  rxn_output = open(outputfilename, 'r')
  number_list = []
  subtraj_list = []
  for line in rxn_output.xreadlines():
    if re.search("<number>",line):
      # pull out the text within the center of the tag
      number_list.append(int(line.strip().split()[1]))
    elif re.search("<subtrajectory>",line):
      # pull out the text within the center of the tag
      subtraj_list.append(int(line.strip().split()[1]))

  rxn_output.close()
  numlines = len(number_list)

  for f in range(numlines):
    # first get the indeces of the trajectory number and subtrajectory
    if counter > max_structures:
      quitting = True
      break
    rxn_number = number_list[f]
    rxn_subtraj = subtraj_list[f]
    # we need to run process_trajectories to pull out 
    # all the trajectory information
    stem = os.path.join(workdir,"proc_traj%d_%d" % (i, f))
    xmltrajfilename = "%s.xml" % (stem,)
    cmd = "process_trajectories -traj %s -index %s -n %d -sn %d \
        -nstride 1 > %s" % (traj_filename, trajindex_filename, 
                            rxn_number, rxn_subtraj, xmltrajfilename)
    print "running command:", cmd
    os.system(cmd)
    # read each trajectory, and pull out the very last frame: 
    # the encounter complex
    trajfile = open(xmltrajfilename,'r')
    trajfilelist = trajfile.readlines()
    trajfile.close()
    lastframelist = trajfilelist[:3] + trajfilelist[-9:]
    lastframename = stem + "_last.xml"
    lastframe = open(lastframename, 'w')
    lastframe.writelines(lastframelist) # write the last frame
    lastframe.close()
    # write the last frame as a pqr file
    pqrfile = os.path.join(workdir, "lig%d_%d.pqr" % (i,f))
    cmd = "xyz_trajectory -mol0 %s -mol1 %s -trajf %s -pqr > %s" % (
        empty, pqrxml1, lastframename, pqrfile)
    print "running command:", cmd
    os.system(cmd)
    # the pqr files must be modified
    modify_pqr(pqrfile)
    # write the receptor pqr
    if i == 0 and f == 0:
      pqrfile = '.'.join(os.path.basename('$PQRXML0').split('.')[:-1])
      cmd = "xyz_trajectory -mol0 %s -mol1 %s -trajf %s -pqr > %s.pqr" \
          % (pqrxml0, empty, lastframename, pqrfile)
      print "running command:", cmd
      os.system(cmd)
      # the pqr files must be modified
      modify_pqr(pqrfile+'.pqr')
      cmd = "pqr2xml < %s.pqr> %s.pqrxml; apbs %s.pqr.in > %s.pqr.out" \
          % (pqrfile, pqrfile, pqrfile, pqrfile)
      print "running command:", cmd
      os.system(cmd)

    os.remove(xmltrajfilename)
    counter += 1
"""

make_fhpd_template = """#!/usr/bin/python

'''
bd_fhpd.py
Takes as arguments a glob of pqr files, for which this script will:
1. Create a directory for each one
2. Prepare all necessary BD files
  - convert to pqrxml
  - reaction files
  - input files
3. Prep the simulation files
4. Run the simulations

Necessary arguments:
- input file template
- receptor files (.pqr, .dx, )
- reaction file template
- number of trajectories per starting structure
'''

import os
import sys
import glob
from string import Template
import random

# CONSTANTS
fhpd_dir = 'fhpd'
input_template_filename = "$INPUT_TEMPLATE_FILENAME"
receptor_pqrxml = "$RECEPTOR_PQRXML"
rxns = "$RXNS"
ntraj = "$NTRAJ"
args = $ARGS

# SCRIPT
input_template_file = open(input_template_filename,'r')
input_template_string=''.join(input_template_file.readlines())
input_template_file.close()


#print input_template
input_template = Template(input_template_string)

# 1. Create directory for each pqr file

if not os.path.exists(fhpd_dir):
  os.mkdir(fhpd_dir)

for arg in args: # for each pqr file
  ligname = os.path.basename(arg).split('.')[0]
  dirname = os.path.join(fhpd_dir, ligname)
  if not os.path.exists(dirname):
    os.mkdir(dirname)

  # 2a. make the PQRXML file
  pqrxml = os.path.join(dirname, ligname+'.pqrxml')
  cmd = "pqr2xml < %s > %s" % (arg, pqrxml)
  print "running command:", cmd
  os.system(cmd)

  # 2b.  make bd input file
  recname = receptor_pqrxml.split('.')[0]
  recdx = os.path.join(recname+'.pqr.dx')
  rxnfile = os.path.join('../..',rxns)
  new_input = input_template.substitute(
      REC=recname, RECDX=recdx, RECPQRXML=receptor_pqrxml, LIG=ligname, 
      NTRAJ=ntraj, RXN=rxnfile, RANDOM=int(10000*random.random()))
  # a link for the pqrxml
  if not os.path.exists(os.path.join(dirname, receptor_pqrxml)): 
                        os.link(receptor_pqrxml, 
                            os.path.join(dirname, receptor_pqrxml))
  # a link for the DX file
  if not os.path.exists(os.path.join(dirname, recdx)): 
                        os.link(recdx, os.path.join(dirname, recdx))
  new_input_filename = os.path.join(dirname,'input.xml')
  new_input_file = open(new_input_filename, 'w')
  new_input_file.write(new_input)
  new_input_file.close()

  # 3. prep the simulations
  os.chdir(dirname)
  cmd = 'bd_top input.xml'
  print "running command", cmd
  os.system(cmd)
  cmd = 'nam_simulation %s-%s-simulation.xml' % (recname,ligname)
  print "running command", cmd
  os.system(cmd)


  os.chdir('../..')

"""

fhpd_consolidate_template = """#!/usr/bin/python

'''
fhpd_consolidate.py
descends into the fhpd file tree, recovers all results .xml files, 
and puts all the results into a single results.xml file in the 
current working directory
'''
import os, sys, glob
import xml.etree.cElementTree as ET # for writing xml files

fhpd_dir = "$FHPD_DIR"
lig_dir_glob = "$LIG_DIR_GLOB"
results_name = "$RESULTS_NAME"

def parse_bd_results(bd_results_filename):
  ''' given a BD results file name, will open the file and extract 
  information about state transitions
  '''
  #bd_results_file = open(bd_results_filename, 'r')
  bd_dict = {}
  if os.path.getsize(bd_results_filename) == 0:
    return bd_dict
  try:
    tree = ET.parse(bd_results_filename)
  except SyntaxError:
    return bd_dict
  root = tree.getroot()
  for tag in root:
    if tag.tag == "reactions":
      reactions = tag
      for tag2 in reactions:
        i = 0
        if tag2.tag == "escaped":
          bd_dict['inf'] = int(tag2.text)
        elif tag2.tag == "completed":
          # need to remove the "rxn" from 
          # the beginning of the site string
          site = tag2[0].text[5:]
          n = tag2[1].text
          #name = outer_state[i] + '_' + str(site)
          bd_dict[site] = int(n)
          i += 1

  return bd_dict

def add_dictionaries(dict1, dict2):
  '''
  adds the values numerically within each dictionary
  NOTE: dict1 is updated and returned BY REFERENCE
  '''
  new_dict = dict1
  for key in dict2.keys():
    if key in dict1.keys():
      dict1[key] += dict2[key]
    else:
      dict1[key] = dict2[key]

  return dict1

def make_new_results_file(bd_results_filename, template_filename, bd_dict):
  tree = ET.parse(template_filename)
  root = tree.getroot()
  for tag in root:
    if tag.tag == "reactions":
      reactions = tag
      for tag2 in reactions:
        i = 0
        if tag2.tag == "escaped":
          tag2.text=str(bd_dict['inf'])
        elif tag2.tag == "completed":
          # need to remove the "rxn" from 
          # the beginning of the site string
          site = tag2[0].text[5:]
          tag2[1].text = str(bd_dict[site])
          i += 1
  results_str = ET.tostring(root)
  results_file = open(bd_results_filename,'w')
  results_file.write(results_str)
  results_file.close()
  return

rxn_dict = {}
globlist = glob.glob(os.path.join(fhpd_dir, lig_dir_glob, results_name))
results_filename = ""
for ligfile in globlist:
  # read the results file
  #print "now reading result file:", ligdir
  results_filename = ligfile
  bd_dict = parse_bd_results(results_filename)
  rxn_dict = add_dictionaries(rxn_dict, bd_dict)

# now we need to write a new results file
#print "rxn_dict:", rxn_dict
assert results_filename, "no results files were read. \
    Possibly something is wrong with the glob inside fhpd_consolidate.py???"
make_new_results_file("results.xml", results_filename, rxn_dict)
"""

default_browndye_params = {
  'root':{
    'protein':'false',
    'solvent':{
      'dielectric':'78',
      'debye-length':'7.86',
      'kT':'1',
    },
    'output':'results.xml',
    'start-at-site':'true',
    'trajectory-file':'traj',
    'include-desolvation-forces':'true',
    'n-trajectories':'10000',
    'n-threads':'1',
    'molecule0': {
      'prefix':'prot0',
      'atoms':'prot0.pqrxml',
      'apbs-grids': {
        'grid':'prot0.dx',
      },
      'solute-dielectric':'2.0',
    },
    'molecule1': {
      'prefix':'prot1',
      'atoms':'prot1.pqrxml',
      'all-in-surface':'false',
      'apbs-grids': {
        'grid':'prot1.dx'
      },
      'solute-dielectric':'2.0',
    },
    'time-step-tolerances': {
      'minimum-dx':'0.20',
    },
    'reactions':'rxns.xml',
    'seed':'111113',
    'n-trajectories-per-output':'1000',
    'n-copies':'200',
    'n-bin-copies':'200',
    'n-steps':'1000000',
    'n-steps-per-output':'1000',
    'max-n-steps':'1000000',
    'n-threads':'1',
    #'min-rxn-coord-file':'min_coord',
  },
}
# for each diffusing molecule, creates this block
default_browndye_molecule_block = {
      'prefix':'prot0',
      'atoms':'prot0.pqrxml',
      'apbs-grids': {
        'grid':'prot0.dx',
      }
}


class dict2xml(object):


    def __init__(self, structure):
        if len(structure) == 1:
            self.doc     = Document()
            rootName    = str(structure.keys()[0])
            self.root   = self.doc.createElement(rootName)

            self.doc.appendChild(self.root)
            self.build(self.root, structure[rootName])

    def build(self, father, structure):
        if type(structure) == dict:
            for k in structure:
                tag = self.doc.createElement(k)
                father.appendChild(tag)
                self.build(tag, structure[k])

        elif type(structure) == list:
            grandFather = father.parentNode
            tagName     = father.tagName
            grandFather.removeChild(father)
            for l in structure:
                tag = self.doc.createElement(tagName)
                self.build(tag, l)
                grandFather.appendChild(tag)

        else:
            data    = str(structure)
            tag     = self.doc.createTextNode(data)
            father.appendChild(tag)

    def display(self):
        print self.doc.toprettyxml(indent="  ")

    def text(self):
        return self.doc.toprettyxml(indent="  ")


def create_ghost_atom_in_pqr(pqr,x,y,z):
  """given a pqr structure, will append a ghost atom at the location 
  x,y,z, where x,y,z are float() values.
  """
  atomid = int(pqr.atoms[-1].index) + 1 # get the last atom index
  resid = int(pqr.atoms[-1].resid) + 1
  print "GHOST atom being added: numbered:", atomid
  ghostatom = pdb.Atom(
      record='ATOM', index=atomid, name="GHO", altloc="", resname="GHO", 
      chain="", resid=resid, icode='', x=x, y=y, z=z, charge='0.0', 
      radius='0.0', occupancy='0.0', beta='0.0', element='')
  pqr.atoms.append(ghostatom)
  pqr.num_atoms += 1
  pqr.num_resids += 1
  return atomid

def prettify(elem):
  """return a pretty-printed xml string for the Element"""
  rough_string = ET.tostring(elem, 'utf-8')
  reparsed = minidom.parseString(rough_string)
  return reparsed.toprettyxml(indent="  ")

def make_rxn_criteria(criteria, pqrs):
  """ takes the rxn criteria pairs and makes ghost atoms in each 
  structure, then populates a rxn file

  criteria format:
    = [[(coords of pqr1), (coords of pqr2), distance], next...]

  example:
    = [[(1.0, 2.0, 3.0), (4.0, 5.0, 6.0), 7.0],]
  """
  # create the xml file
  roottag = ET.Element("roottag")
  first_state = ET.SubElement(roottag, "first-state")
  first_state.text = "start"
  reactions = ET.SubElement(roottag, "reactions")
  reaction_list = []
  #pqrs = map(copy.deepcopy, raw_pqrs)
  ghost0_coords_already_used = []
  ghost0_ids = {}
  ghost1_coords_already_used = []
  ghost1_ids = {}
  for i in range(len(criteria)):
    #coord0 = criteria[i][0]
    #coord1 = criteria[i][1]
    radius = criteria[i]['radius']
    # We don't have any duplicate ghost atoms
    # check if there is already a ghost atom in the molecule
    if (criteria[i]['centerx'], criteria[i]['centery'], 
        criteria[i]['centerz']) in ghost0_coords_already_used:
      ghost0_id = ghost0_ids[(criteria[i]['centerx'], criteria[i]['centery'],
                              criteria[i]['centerz'])]
    else: # then the ghost atom doesn't already exist
      # the star spreads the coords out into individual arguments
      ghost0_id = create_ghost_atom_in_pqr(
          pqrs[0], criteria[i]['centerx'], criteria[i]['centery'], 
          criteria[i]['centerz']) 
      ghost0_coords_already_used.append((
          criteria[i]['centerx'], criteria[i]['centery'], 
          criteria[i]['centerz']))
      ghost0_ids[(
          criteria[i]['centerx'], criteria[i]['centery'], 
          criteria[i]['centerz'])] = ghost0_id

    # check the ligand for already existing ghost atoms
    # check if there is already a ghost atom in the molecule
    if (criteria[i]['ligx'], criteria[i]['ligy'], criteria[i]['ligz']) \
        in ghost1_coords_already_used:
      ghost1_id = ghost1_ids[(criteria[i]['ligx'], criteria[i]['ligy'], 
                              criteria[i]['ligz'])]
    else: # then the ghost atom doesn't already exist
      ghost1_id = create_ghost_atom_in_pqr(
          pqrs[1], criteria[i]['ligx'], 
          criteria[i]['ligy'], criteria[i]['ligz'])
      ghost1_coords_already_used.append((
          criteria[i]['ligx'], criteria[i]['ligy'], 
          criteria[i]['ligz']))
      ghost1_ids[(
          criteria[i]['ligx'], criteria[i]['ligy'], 
          criteria[i]['ligz'])] = ghost1_id


    siteid = criteria[i]['siteid']
    index = criteria[i]['index']
    # now add to the xml file
    reaction_list.append(ET.SubElement(reactions, "reaction"))
    name=ET.SubElement(reaction_list[-1],"name")
    name.text = "%s_%s" % (siteid,index)
    state_before = ET.SubElement(reaction_list[-1], "state-before")
    state_before.text = "start"
    state_after = ET.SubElement(reaction_list[-1], "state-after")
    state_after.text = "end"
    criterion = ET.SubElement(reaction_list[-1], "criterion")
    n_needed = ET.SubElement(criterion, "n-needed")
    n_needed.text = "1"
    pair = ET.SubElement(criterion, "pair")
    atoms = ET.SubElement(pair, "atoms")
    atoms.text = "%d %d" % (ghost0_id, ghost1_id)
    distance = ET.SubElement(pair, "distance")
    distance.text = "%f" % (radius,)

  criteria_xml = prettify(roottag)
  return criteria_xml

def pqr2xml(pqrfile, pqr2xml_program='pqr2xml'):
  """simply runs the pqr2xml program that comes with Browndye."""
  #print "pqr2xml_program:", pqr2xml_program
  # get everything but the extension
  no_ext = os.path.splitext(pqrfile)[0]
  xmlfile = no_ext + '.pqrxml'
  command = '%s < %s > %s' % (pqr2xml_program, pqrfile, xmlfile)
  if verbose: print "now running command:", command
  result = os.system(command) # run the pqr2xml program
  if result != 0: raise Exception, "There was a problem running pqr2xml"
  return xmlfile


def write_browndye_input(pqrs,seekrcalc,criteria,work_dir='.',browndye_bin='', 
                         start_at_site='true',fhpd_mode=False):
  """generates a Browndye input file"""

  counter = 0
  input_xml = default_browndye_params
  debyes = []
  #apbs_settings = seekrcalc.browndye.apbs #settings['apbs_settings']
  #inputgen_settings = seekrcalc.browndye.apbs.inputgen #settings[
  #   'inputgen_settings']
  rxn_criteria = make_rxn_criteria(criteria,pqrs)
  # modifies the molecules to contain ghost atoms 
  # (for BrownDye rxn criteria)
  pqrxmls = []
  for pqr in pqrs: # for each molecule in pqr format
    prefix = pqr.struct_id # name of the molecule
    #print "PREFIX", prefix
    pqrfile = os.path.join(work_dir, prefix+'.pqr')
    pqr.save(pqrfile,pqr=True,endmdl=False)
    #print "pqrfile:", pqrfile
    # get the electrostatic grid and debye length for the molecule
    dxfile, debye = apbs.main(seekrcalc, pqrfile, fhpd_mode=fhpd_mode)
    debyes.append(debye)
    # call the pqrxml program using the Browndye software suite
    pqrxmlfile = pqr2xml(pqrfile, pqr2xml_program=os.path.join(browndye_bin, 
                                                               'pqr2xml'))
    pqrxmls.append(pqrxmlfile)
    # create a copy of the molecule block, 
    # keep the default solute dielectric
    molecule_xml = deepcopy(default_browndye_molecule_block)
    molecule_xml['apbs-grids']['grid']=os.path.basename(dxfile)
    molecule_xml['prefix'] = prefix
    molecule_xml['atoms'] = os.path.basename(pqrxmlfile)
    molecule_xml['solute-dielectric'] = '2.0'
    # copy the molecule tree over to the input xml
    input_xml['root']['molecule%d'%counter] = molecule_xml

    counter += 1

  #for pqr in pqrs:
    #pqr.save() # save the pqr files containing the ghost molecules
  rxn_file = open(os.path.join(work_dir,RXN_FILENAME), 'w')
  rxn_file.write(rxn_criteria)
  rxn_file.close()
  # give each one a random seed to make the simulation rounds different
  input_xml['root']['seed'] = int(random.random() * 10000)
  input_xml['root']['protein'] = seekrcalc.browndye.ligand_is_protein
  # set the debye length
  input_xml['root']['solvent']['debye-length'] = str(max(map(int,map(float,
                                                                     debyes))))
  input_xml['root']['reactions'] = RXN_FILENAME
  #float(settings['temperature']) / DEFAULT_TEMP
  input_xml['root']['solvent']['kT'] = seekrcalc.master_temperature 
  input_xml['root']['n-threads'] = seekrcalc.browndye.num_threads
  input_xml['root']['n-trajectories'] = seekrcalc.browndye.prods_per_anchor
  input_xml['root']['start-at-site'] = start_at_site
  if fhpd_mode:
    input_xml['root']['seed'] = "$RANDOM"
    input_xml['root']['include-desolvation-forces'] = "true"
    # remove the apbs-grid because we have no desolvation forces
    del input_xml['root']['molecule1']['apbs-grids']
    input_xml['root']['molecule0']['prefix'] = "$REC"
    input_xml['root']['molecule0']['atoms'] = "$REC.pqrxml"
    input_xml['root']['molecule0']['apbs-grids']['grid'] = "$RECDX"
    input_xml['root']['molecule1']['prefix'] = "$LIG"
    input_xml['root']['molecule1']['atoms'] = "$LIG.pqrxml"
    input_xml['root']['reactions'] = "$RXN"
    input_xml['root']['n-trajectories'] = "$NTRAJ"
  
  # generate xml text for the file
  input_text = dict2xml(input_xml).text()
  input_xml = {}
  input_file = open(os.path.join(work_dir,INPUT_FILENAME), 'w')
  input_file.write(input_text) # write an xml file for the input to bd
  input_file.close()

  return pqrxmls

def shape_intersects_sphere_stupid(sphere, shape):
  """ Determines whether a shape intersects the shell of a sphere.
  It uses a stupid method to do this which isn't perfect but should be 
  good enough for this purpose (because the sphere is big compared to 
  the box) Stupid method: iterates through all the corner points. If 
  they are all inside or all outside the sphere, the function returns 
  False, otherwise True

  sphere: a dictionary containing 'x','y','z', & 'radius' indeces
  shape: a list of 3-tuples of the points defining a shape
  """
  consensus = None
  for point in shape:
    within = True
    # see if the distance between the point and sphere center is greater
    # than the sphere radius
    if math.sqrt((sphere['x'] - point[0])**2 + (sphere['y'] - point[1])**2 
                  + (sphere['z'] - point[2])**2) > sphere['radius']:
      within = False

    if consensus not in [within, None]:
      # we do not have consensus, therefore some points are within 
      # while others are without, therefore shapes do intersect
      return True
    else:
      consensus = within
  return False # there is consensus, therefore no intersection

def make_empty_pqrxml(filename):
  """Create the empty.pqrxml file."""
  emptyxml = open(filename, 'w')
  emptyxml.write(empty_pqrxml_template)
  emptyxml.close()

def build_bd(seekrcalc):
  """build all structures and necessary files for BD calculations

  takes a list of pqr filenames, reaction coordinate pairs/distances
  """

  if verbose: 
    print '\n', '#'*40, "\n \Now creating BD files using bd.py\n", '#'*40

  parser = pdb.Big_PDBParser()
  rec_struct = parser.get_structure(
      'bd_receptor_dry_pqr', seekrcalc.browndye.rec_dry_pqr_filename, pqr=True)

  #milestone_pos_rot_list = settings['milestone_pos_rot_list'] 
  # # NOTE: may want to clean up code referring to this variable
  
  """ # TODO: marked for removal because feature probably not needed
  if settings['starting_conditions'] == 'configs':
    lig_configs = settings['lig_configs']
  elif settings['starting_conditions'] == 'spheres':
    # then generate the anchors from 
    # random positions/orientations in a sphere
    # NOTE: there is a random orientation 
    # function in positions_orient.py
    pass
  else:
    raise Exception, "option not allowed: %s" % settings['starting_conditions']
  """
  
  #bd_file_paths = settings['bd_file_paths']
  #browndye_bin = settings['browndye_bin_dir']
  #if not os.path.exists(empty_pqrxml)
  
  bd_configs = []

  # fill the b_surface folder
  # get the first config of the ligand
  lig_config = seekrcalc.browndye.starting_lig_config
  lig_center = pdb.center_of_mass(lig_config)
  pqrs = [copy.deepcopy(rec_struct), copy.deepcopy(lig_config)]
  pqrs[1].struct_id='bd_ligand'
  

  #for site in settings['b_surface_ending_surfaces']:
  b_surface_criteria = []
  starting_surfaces = []
  for milestone in seekrcalc.milestones:
    if milestone.end:
      # add every site to the criteria list
      b_surface_criteria.append({
          'centerx':milestone.center_vec[0], 'centery':milestone.center_vec[1],
          'centerz':milestone.center_vec[2], 'ligx':lig_center[0], 
          'ligy':lig_center[1],'ligz':lig_center[2], 'radius':milestone.radius, 
          'index':milestone.index, 'siteid':milestone.siteid})
    if milestone.bd:
      starting_surfaces.append({
          'site':milestone.siteid, 'radius':milestone.radius, 
          'index':milestone.index})
  
  print "bsurface_criteria:", b_surface_criteria
  # write input for this part
  b_surface_pqrxmls = write_browndye_input(
      pqrs, seekrcalc, b_surface_criteria, 
      work_dir=seekrcalc.browndye.b_surface_path, 
      browndye_bin=seekrcalc.browndye.browndye_bin, 
      start_at_site='false', fhpd_mode = False)
  
  """
  for bd_file_path in bd_file_paths:
    if not bd_file_path:
      pass
      #bd_configs.append(None)
      #continue # then don't do BD for this portion
    #print "bd_file_path:", bd_file_path
    # reading the file tree to get the index of every existing folder
    anchor_folder_name = bd_file_path.split('/')[-2]
    #print "anchor_folder_name:", anchor_folder_name
    # getting the index out of the folder name
    folder_index = anchor_folder_name.split('_')[1]
    bd_configs.append(int(folder_index))
  """
  counter = 0 # the index of the loop itself
  
    
  for milestone in seekrcalc.milestones:
    if milestone.bd:
      lig_config = milestone.config
      bd_file_path = os.path.join(seekrcalc.project.rootdir, 
                                  milestone.directory, 'bd')
      print "bd_file_path:", bd_file_path
      pqrs = [copy.deepcopy(rec_struct), copy.deepcopy(lig_config)]
      pqrs[1].struct_id='bd_ligand'
      bd_needed = True
      
      """
      m=0
      # m will represent the proper milestone index
      for m in range(len(milestone_pos_rot_list)):
        if milestone_pos_rot_list[m][0].index == i:
          break
      """
      criteria = []
      for surface in starting_surfaces:
        if surface['site'] == milestone.siteid:
          proper_radius = milestone.bd_adjacent.radius
          proper_index = milestone.bd_adjacent.index
        else:
          proper_radius = surface['radius']
          proper_index = surface['index']
        criteria.append({
            'centerx':milestone.center_vec[0], 
            'centery':milestone.center_vec[1], 
            'centerz':milestone.center_vec[2], 
            'ligx':lig_center[0], 'ligy':lig_center[1], 
            'ligz':lig_center[2], 'radius':proper_radius, 
            'index':proper_index, 'siteid':milestone.siteid})
      """
      #[[(31.121, 37.153, 35.253), (38.742, 51.710, 68.137), 9.0],] 
      # a list of all reaction criteria
      criteria = []
      for site in settings['starting_surfaces']:
        #site_center = [site['x'], site['y'],site['z']]
        # NOTE: at this time, only spherical reaction 
        # criteria are allowed in BrownDye
        #radius = site['radius']
        # then its the same site, we need to go one shell in
        if milestone_pos_rot_list[m][0].siteid == site['siteid']:
          # the radius in the same site
          proper_radius = site['inner_radius']
          proper_index = site['inner_index']
        # this is a different site, choose the same radius as starting
        else:
          proper_radius = site['outer_radius']
          proper_index = site['outer_index']
        # add every site to the criteria list
        criteria.append({
            'centerx':site['x'], 'centery':site['y'], 
            'centerz':site['z'], 'ligx':lig_center[0], 
            'ligy':lig_center[1], 'ligz':lig_center[2], 
            'radius':proper_radius, 'index':proper_index, 
            'siteid':site['siteid']})
      """
      #print "pqrs:", pqrs, 'criteria:', criteria
    
    
      site_pqrxmls = write_browndye_input(
          pqrs, seekrcalc, criteria, work_dir=bd_file_path, 
          browndye_bin=seekrcalc.browndye.browndye_bin,fhpd_mode=True)
      # make BD preparation scripts extract_bd_frames.py and bd_fhpd.pyp
      # Write the script to extract all frames from the 
      # successful b_surface bd trajectories
    
      extract_bd_frames_dict = {
          'TRAJDIR':"../../b_surface", 'WORKDIR':"./trajs", 
          'PQRXML0':os.path.basename(b_surface_pqrxmls[0]), 
          'PQRXML1':os.path.basename(b_surface_pqrxmls[1]), 
          'EMPTY':empty_pqrxml, 
          'SITENAME':'%s_%s' % (milestone.siteid, milestone.index), 
          'NUMBER_OF_TRAJS':seekrcalc.browndye.num_threads}
      extract_bd_frames = Adv_template(extract_bd_frames_template, 
                                       extract_bd_frames_dict)
      extract_file = open(os.path.join(bd_file_path,
                                      "extract_bd_frames.py"), 'w')
      # write an xml file for the input to bd
      extract_file.write(extract_bd_frames.get_output())
      extract_file.close()
      # construct the FHPD distribution prep scripts
      make_fhpd_dict = {
          'INPUT_TEMPLATE_FILENAME':'input.xml', 
          'RECEPTOR_PQRXML':os.path.basename(b_surface_pqrxmls[0]), 
          'RXNS':'rxns.xml', 'NTRAJ':seekrcalc.browndye.fhpd_numtraj, 
          'ARGS':"glob.glob(os.path.join('./trajs','lig*.pqr'))"}
      # NOTE: should change NTRAJ to be consistent with the number of 
      # reaction events in the b_surface phase
      make_fhpd = Adv_template(make_fhpd_template,make_fhpd_dict)
      make_fhpd_file = open(os.path.join(bd_file_path,"make_fhpd.py"), 'w')
      # write an xml file for the input to bd
      make_fhpd_file.write(make_fhpd.get_output()) 
      make_fhpd_file.close()
      # Consolidate all result files from the FHPD simulations into 
      # one large results.xml file
      fhpd_consolidate_dict = {'FHPD_DIR':"fhpd", 'LIG_DIR_GLOB':"lig*/", 
                               'RESULTS_NAME':'results.xml'}
      fhpd_consolidate = Adv_template(fhpd_consolidate_template,
                                      fhpd_consolidate_dict)
      fhpd_consolidate_file = open(os.path.join(bd_file_path,
                                   "fhpd_consolidate.py"), 'w')
      # write an xml file for the input to bd                             
      fhpd_consolidate_file.write(fhpd_consolidate.get_output()) 
      fhpd_consolidate_file.close()
      make_empty_pqrxml(os.path.join(bd_file_path, 'empty.pqrxml'))

    counter += 1




class Test_bd_functions(unittest.TestCase):
  # several test cases to ensure the functions 
  # in this module are working properly
  def test_main(self): # test this function
    pass

  def test_write_browndye_input(self):
    print "Warning: unit tests for function 'write_browndye_input' \
          have not been properly implemented."
    # It may to be too difficult to run this test 
    # as it needs complicated objects
    return 

  def test_dict_to_xml(self):
    test_dict = {
        'root':{'layer1a':{'layer2a':'setting1',},
        'layer1b':'setting2','layer1c':'setting3'}}
    test_xml = '<?xml version="1.0" ?>\n<root>\n  \
        <layer1c>setting3</layer1c>\n  <layer1b>setting2</layer1b>\n  \
        <layer1a>\n    <layer2a>setting1</layer2a>\n  </layer1a>\n</root>\n'
    self.assertEqual(dict2xml(test_dict).text(),test_xml)

  def test_pqr2xml(self):
    
    browndye_bin = os.environ['BROWNDYE_BIN']
    self.XML = pqr2xml(test_pqr_filename, 
                       pqr2xml_program=os.path.join(browndye_bin, 'pqr2xml'))
    fileexists = os.path.exists(self.XML)
    self.assertTrue(fileexists)

  def test_make_rxn_criteria(self):

    test_xml = """<?xml version="1.0" ?>\n<roottag>\n  \
        <first-state>start</first-state>\n  <reactions>\n    \
        <reaction>\n      <name>0_0</name>\n      \
        <state-before>start</state-before>\n      \
        <state-after>end</state-after>\n      <criterion>\n        \
        <n-needed>1</n-needed>\n        <pair>\n          \
        <atoms>4323 4324</atoms>\n          \
        <distance>7.000000</distance>\n        </pair>\n      \
        </criterion>\n    </reaction>\n  </reactions>\n</roottag>\n"""
    self.assertEqual(test_xml, make_rxn_criteria(test_rxns, 
                                                 [test_pqr,test_pqr]))

  #def test_anchor_touching_starting_surface(self):
  #  pass # I don't think this function is used anywhere


  def test_shape_intersects_sphere_stupid(self):
    sphere = {'x':0.0,'y':0.0,'z':0.0,'radius':10.0}
    shape_in = [(1.0,1.0,1.0),(-1.0,1.0,1.0),(1.0,-1.0,1.0),(1.0,1.0,-1.0)]
    shape_out = [(101.0,101.0,101.0),(99.0,101.0,101.0),
                 (101.0,99.0,101.0),(101.0,101.0,99.0)]
    shape_intersect = [(1.0,1.0,1.0),(10.0,10.0,10.0),
                      (10.0,1.0,1.0),(1.0,10.0,1.0)]
    self.assertEqual(False, 
                     shape_intersects_sphere_stupid(sphere,shape_in))
    self.assertEqual(False, 
                     shape_intersects_sphere_stupid(sphere,shape_out))
    self.assertEqual(True, 
                     shape_intersects_sphere_stupid(sphere,shape_intersect))


if __name__=='__main__':
  if __name__ == "__main__" and 'INPUTGEN' not in os.environ:
    print "In order to run unit tests in bd.py, please set the BROWNDYE_BIN \
          environmental variable to the directory containing the BrownDye \
          programs"
    exit()
  print "Running unit tests for bd.py"
  test_holo = parser.get_structure('small','../test/test_tiny.pdb')
  test_inputgen_filename = os.environ['INPUTGEN']
  test_pqr_filename = "../test/1cbj.pqr"
  #old_test_rxns = [[(1.0, 2.0, 3.0), (4.0, 5.0, 6.0), 7.0],]
  test_rxns = [{'centerx':1.0, 'centery':2.0, 'centerz':3.0, 
                'ligx':4.0, 'ligy':5.0, 'ligz':6.0, 
                'radius':7.0, 'index':0, 'siteid':0}]
  test_pqr = parser.get_structure('pqr', test_pqr_filename, pqr=True)
  unittest.main() # run tests of all functions


