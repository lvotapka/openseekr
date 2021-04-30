#!/usr/bin/python

"""
apbs.py

creates the necessary files to run an electrostatic simulation 
using APBS

"""
import os
import sys
from math import sqrt

import shutil #, math, subprocess #, make_fxd
import pdb2 as pdb
import unittest
import re

from adv_template import *


verbose = True

parser = pdb.Big_PDBParser()

# get the path to this script
self_path = os.path.dirname(os.path.realpath(__file__)) 

apbs_input_template_location = os.path.join(self_path, 
                                            "apbs_input.template")

if __name__ == "__main__" and "INPUTGEN" not in os.environ:
  print "In order to run unit tests in apbs.py, \
      please set the INPUTGEN environmental variable"
  exit()

if __name__ == "__main__":
  test_inputgen_location = os.environ["INPUTGEN"]
  test_apbs_location = "apbs"
  test_pqr_filename = "../test/1cbj.pqr"

k_B = 1.380649e-23 # m^2 kg per s^2 per K
coulombs_per_proton = 1.602177e-19 # C
permittivity = 8.854188e-12 * 78 # C^2 per N per m^2
molecules_per_m3_per_moles_per_liter = 6.022e26
angstroms_per_m = 1e10

default_apbs_params = {
"pqr":"",
"dimx":"65",
"dimy":"65",
"dimz":"65",
"cglenx":"100.0000",
"cgleny":"100.0000",
"cglenz":"100.0000",
"fglenx":"65.0000",
'fgleny':'65.0000',
"fglenz":"65.0000",
"boundary_condition":"sdh",
"lpbe_npbe":"lpbe",
"solute_dielec":"2.0",
"solvent_dielec":"78.5400",
"temp":"310.0",
"stem":"pot",
}

default_inputgen_settings = {
  "fadd":"60",
  "gmemceil":"64000",
  "resolution":"0.5",
  "ionic_str":"0.15",
  "cfac": "4.0",
}


def make_apbs_input_using_inputgen(
    inputgen_filename, pqr_filename, fadd=60, 
    cfac=4.0, gmemceil=64000, resolution=0.5, ionic_str=0.15):
  """makes an apbs input file given a pqr file & other parameters. 
  See Inputgen.py in PDB2PQR documentation for descriptions 
  of other arguments
  """
  pqr_basename= os.path.basename(pqr_filename)
  pqr_abspath=os.path.dirname(os.path.abspath(pqr_filename))
  pre_ext = (pqr_basename.split("."))[0]
  runstring = "python %s --potdx --fadd=%s --cfac=%s --space=%s \
              --gmemceil=%s --istrng=%s %s" % (inputgen_filename, fadd, cfac, 
                                               resolution, gmemceil, ionic_str, 
                                               pqr_basename,)
  olddir = os.curdir
  os.chdir(os.path.dirname(pqr_filename))
  print "Now creating APBS input file using command:", runstring
  os.system(runstring)
  os.chdir(olddir)
  input_filename = os.path.join(pqr_abspath+"/"+ pre_ext+".in")
  return input_filename

def scrape_inputfile(input_filename):
  """NOTE: only takes out the dime, pdime, cglen, 
  and fglen parameters from an APBS input file.
  """
  dimestring = pdimestring = cglenstring = fglenstring = None
  infile = open(input_filename,"r")
  for line in infile:
    if re.search(" dime", line) and not dimestring:
      dimestring = line
    if re.search("pdime", line) and not pdimestring:
      pdimestring = line
    if re.search("cglen", line) and not cglenstring:
      cglenstring = line
    if re.search("fglen", line) and not fglenstring:
      fglenstring = line
  infile.close()
  if pdimestring: 
    raise Exception, "Parallel-run dx files not yet implemented..."
  apbs_params = {} # a dictionary containing what we scraped outta here
  dime_list = dimestring.strip().split()
  cglen_list = cglenstring.strip().split()
  fglen_list = fglenstring.strip().split()
  apbs_params["dimx"],apbs_params["dimy"],apbs_params["dimz"] = dime_list[1:]
  apbs_params["cglenx"],apbs_params["cgleny"],apbs_params["cglenz"] \
  = cglen_list[1:]
  apbs_params["fglenx"],apbs_params["fgleny"],apbs_params["fglenz"] \
  = fglen_list[1:]
  return apbs_params

def make_apbs_input_using_template (new_apbs_params, 
                                    apbs_file_location="apbs.in"):
  # create empty directory
  apbs_params = {} 
  # populate with default values
  apbs_params.update(default_apbs_params) 
  # populate with new parameters
  apbs_params.update(new_apbs_params) 
  # fill parameters into the template to make apbs file
  apbs_input = File_template(apbs_input_template_location, apbs_params)
  # save the apbs input file
  apbs_input.save(apbs_file_location) 
  return

def run_apbs (apbs_filename, input_filename, 
              pqr_filename, std_out="apbs.out"):
  """runs apbs using a given input file "input_filename" and writes all 
  standard output to 'std_out'.
  """
  rundir = os.path.dirname(input_filename)
  print "copying file: %s to directory: %s" % (
      pqr_filename, os.path.join(rundir,os.path.basename(pqr_filename)))
  if os.path.abspath(os.path.dirname(pqr_filename)) != os.path.abspath(rundir):
    shutil.copyfile(pqr_filename, 
                    os.path.join(rundir,os.path.basename(pqr_filename)))
  pqr_filename = os.path.basename(pqr_filename)
  input_filename = os.path.basename(input_filename)
  std_out = os.path.basename(std_out)
  # string to run apbs
  runstring = "%s %s > %s" % (apbs_filename, input_filename, std_out) 
  
  curdir = os.getcwd()
  # inputgen will automatically make this .dx file
  dx_filename = pqr_filename + ".dx"
  if os.path.exists(dx_filename):
    print "DX file found:", dx_filename
    print "skipping APBS run..."
  else:
    # we want to run APBS in the directory
    os.chdir(rundir) 
    print "running command:", runstring
    # execute the string
    result = os.system(runstring) 
    # then an error occured
    if result != 0: raise Exception, "There was a problem running APBS" 
    os.chdir(curdir)
  
  return dx_filename # return the name of the dx file

def get_debye_length(seekrcalc, apbs_std_outfilename):
  """Will parse an apbs stdout file to look for the Debye length."""
  debye_string = re.compile("Debye length")
  debye_list = [] # a list of numbers that will be returned
  for line in open(apbs_std_outfilename, "r"):
    m = re.search(debye_string, line)
    if m: # then we've found a line
      number_obj = re.search("[0-9.]+", line).group()
      if number_obj == "0": 
        print "ALERT: Debye length of zero found. This may mean that your PQR \
            file has a net charge that is NOT zero, or that your ion \
            concentration was zero..."
      else:
        debye_list.append(number_obj)
  if len(debye_list) < 1:
    print "APBS was unable to find a Debye length. \
        Now computing using a formula..."
    denominator = 0.0
    for ion in seekrcalc.browndye.apbs.ions:
      denominator += ion.concentration * molecules_per_m3_per_moles_per_liter \
      * (ion.charge*coulombs_per_proton)**2
    inside = float(
        permittivity * k_B * seekrcalc.master_temperature / denominator)
    debye_length = sqrt(inside) * angstroms_per_m
    debye_list.append(debye_length)  
    print "A Debye length of", debye_length, "was found."
  assert len(debye_list) > 0, "Debye length not found in APBS output: %s. \
    Please ensure that APBS calculation was completed properly and that the \
    correct output file was specified."
  # take the first member of it by default. There may be a better way 
  # for this but all outputs seem to be the same
  return debye_list[0]

def flatten_ion_list(apbs_settings):
  if "ions" not in apbs_settings.keys(): return apbs_settings
  ion_list = apbs_settings.pop("ions")
  for ion in ion_list:
    key = ion["key"]
    apbs_settings["%sconc" % key] = ion["concentration"]
    apbs_settings["%scrg" % key] = ion["charge"]
    apbs_settings["%srad" % key] = ion["radius"]
  return apbs_settings

def main(seekrcalc, pqr_filename,fhpd_mode=False):
  #user_settings = {}
  #user_settings.update(default_inputgen_settings)
  #user_settings.update(inputgen_settings)
  #apbs_settings = flatten_ion_list(apbs_settings)
  #if apbs_settings['ion1conc']: 
  #   user_settings['ionic_str'] = apbs_settings['ion1conc']
  inputgen = seekrcalc.browndye.apbs.inputgen
  apbs = seekrcalc.browndye.apbs
  if apbs.ions[0].concentration and not inputgen.ionic_str:
    inputgen.ionic_str = apbs.ions[0].concentration
  # make APBS input file using inputgen (enabled)
  inputgen_location = inputgen.executable
  apbs_location = apbs.executable
  if not os.path.exists(inputgen_location):
    print "provided inputgen.py location does not exist. Please modify \
        browndye.inputgen.executable value."
    exit()
  input_filename = make_apbs_input_using_inputgen(
      inputgen_location, pqr_filename, fadd=inputgen.fadd, cfac=inputgen.cfac, 
      gmemceil=inputgen.gmemceil, resolution=inputgen.resolution, 
      ionic_str=inputgen.ionic_str)
  # make APBS input file using template (disabled)
  # make DX grids
  apbs_out=pqr_filename+".out" # make a default apbs output file
  # use the inputgen-generated file to 
  # make our own, more customized file
  apbs_params = scrape_inputfile(input_filename)
  #if fhpd_mode:
  apbs_params["pqr"] = apbs_params["stem"] = os.path.basename(pqr_filename)
  #else:
    #apbs_params['pqr'] = apbs_params['stem'] = pqr_filename
  for i, ions in enumerate(apbs.ions):
    print "ion:", ions.concentration, ions.charge, ions.radius
    apbs_params["ion%dconc" % i] = ions.concentration
    apbs_params["ion%dcrg" % i]  = ions.charge
    apbs_params["ion%drad" % i]  = ions.radius
  
  new_input_filename = pqr_filename + ".in"
  dx = pqr_filename + ".dx"
  print "apbs_params:", apbs_params
  make_apbs_input_using_template(apbs_params, new_input_filename)
  #if not fhpd_mode:
  # save the electrostatic grid
  run_apbs(apbs_location, new_input_filename, pqr_filename, std_out=apbs_out)
  # find the Debye length
  debye = get_debye_length(seekrcalc, apbs_out)
  return dx, debye

def is_number(s):
  """returns True if the string 's' can be converted to a float/int, 
  False otherwise.
  """
  try:
    float(s)
    return True
  except ValueError:
    return False

class Test_apbs_functions(unittest.TestCase):
  # several test cases to ensure the functions 
  # in this module are working properly
  # test whether the apbs input file has been created properly
  def test_make_apbs_input_using_inputgen(self): 
    #print 'test pqr filename', test_pqr_filename
    self.APBS_inp = make_apbs_input_using_inputgen(
        test_inputgen_location, test_pqr_filename) # get file location
    #print "self.APBS_inp", self.APBS_inp
    fileexists = os.path.exists(self.APBS_inp)
    self.assertTrue(fileexists)

  # test whether the apbs input file has been created properly
  def test_make_apbs_input_using_template(self):
    make_apbs_input_using_template({},"/tmp/test_apbs.in")
    fileexists = os.path.exists("/tmp/test_apbs.in")
    # if it exists, then that's good
    self.assertTrue(fileexists) 

  # test whether apbs is running properly
  def test_run_apbs(self): 
    # get file location
    self.APBS_inp = make_apbs_input_using_inputgen(
        test_inputgen_location, os.path.abspath(test_pqr_filename)) 
    self.APBS_inp2 = "/tmp/input2.in"
    self.inp_dict = scrape_inputfile(self.APBS_inp)
    self.inp_dict["pqr"] = self.inp_dict["stem"] \
    = os.path.abspath(test_pqr_filename)
    make_apbs_input_using_template(self.inp_dict, self.APBS_inp2)
    run_apbs(test_apbs_location, self.APBS_inp2, test_pqr_filename)
    self.APBS_dx = test_pqr_filename + ".dx"
    fileexists = os.path.exists(self.APBS_dx)
    self.assertTrue(fileexists)

  def test_is_number(self): # test the is_number function
    self.assertTrue(is_number("0"))
    self.assertTrue(is_number("3.14"))
    self.assertTrue(is_number("2.0e-8"))
    self.assertFalse(is_number("foobar"))

  def test_get_debye_length(self):
    testfile1 = open("/tmp/debye_test1","w") # file with numbers
    testfile1.writelines(["CALCULATION #1: MULTIGRID\n",
				"Setting up problem...\n",
  				"Vpbe_ctor:  Using max ion radius (0 A) for exclusion function\n",
  				"Debye length:  99.23 A\n",
  				"Current memory usage:  731.506 MB total, 731.506 MB high water\n",
  				"Using cubic spline charge discretization.\n",])
    testfile1.close()
    testfile2 = open("/tmp/debye_test2","w") # file with nothing
    testfile2.close()
    # should return the numeric value 99.32
    result1 = get_debye_length("/tmp/debye_test1") 
    self.assertEqual(result1, "99.23")
    # this is an empty file, so the function should throw an error
    self.assertRaises(AssertionError, get_debye_length, "/tmp/debye_test2") 

  def test_scrape_inputfile(self):
    testfile1 = open("/tmp/scrape_test1","w") # file with numbers
    testfile1.writelines(["""    dime 129 129 193
    cglen 80.2842 77.5999 116.9345
    fglen 67.2260 65.6470 88.7850
    """,])
    testfile1.close()
    test_params = {"dimx":"129", "dimy":"129", "dimz":"193","cglenx":"80.2842",
                   "cgleny":"77.5999", "cglenz":"116.9345", "fglenx":"67.2260", 
                   "fgleny":"65.6470", "fglenz":"88.7850"}
    self.assertEqual(test_params, scrape_inputfile("/tmp/scrape_test1"))


if __name__=="__main__":
  print "Running unit tests for apbs.py"
  unittest.main() # run tests of all functions
