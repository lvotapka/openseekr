'''
Functions and objects intended to run the reversal and forward stages of 
a SEEKR calculation.

Created on June 15, 2018

@author: lvotapka
'''

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import sys, os, time
from sys import stdout
from seekr import amber
from copy import deepcopy

from seekrplugin import SeekrForce
from adv_template import Adv_template

import cPickle as pickle

verbose = True

MAX_REVERSE_ITER = 100000 # 200 ns
MAX_FORWARD_ITER = 100000 # 200 ns

def make_box_info(box_vectors):
  '''Converts a 3x3 matrix of box vectors to x,y,z,alpha,beta,gamma format for 
  CPPTRAJ input.
  Input:
   - box_vectors: A 3x3 matrix representing the triclinic box vectors
  Output:
   - box_info_string: a string formatted for CPPTRAJ 'box' command
  '''
  box_length = box_vectors[0][0].value_in_unit(angstroms) # extract the number without units
  box_info_string = 'x %f y %f z %f alpha 109.4712190 beta 109.4712190 gamma 109.4712190' % (box_length, box_length, box_length)
  return box_info_string

def autoimage_traj(parm_name, trajin_name, trajout_name, box_info, cpptraj_script_location, cpptraj_exe='cpptraj'):
  '''Runs the CPPTRAJ autoimage command for a triclinic box simulation.
  Input:
   - parm_name: a string representing the filename of a .prmtop or .parm7 AMBER
       parameter/topology file
   - trajin_name: a non-imaged trajectory to load for imaging
   - box_info: a string representing the triclinic box in x,y,z,alpha,beta,gamma
       format
   - cpptraj_script_location: string for the location to write the cpptraj script
   - cpptraj_exe: an optional string representing the OS command to run CPPTRAJ
  Output:
   - None
   '''
  cpptraj_template = '''parm $PARMFILE
box $BOX_INFO
trajin $TRAJIN
autoimage
trajout $TRAJOUT
go
quit
'''
  cpptraj_dict = {'PARMFILE':parm_name, 'TRAJIN':trajin_name, 'TRAJOUT':trajout_name, 'BOX_INFO':box_info} # define template dictionary
  cpptraj_script = Adv_template(cpptraj_template, cpptraj_dict) # fill in the values into the template from the dictionary
  extract_file = open(cpptraj_script_location, 'w') # open the script for writing
  extract_file.write(cpptraj_script.get_output()) # write a cpptraj script
  extract_file.close()
  os.system("%s < %s" % (cpptraj_exe, cpptraj_script_location)) # run cpptraj
  return

def create_spherical_seekr_force(seekrcalc, milestone, system, end_on_middle_crossing):
  '''create the SEEKR 'force' even though it's more of a monitor than a force.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to run the simulation for
   - system: the OpenMM system to add the SEEKR 'force' to
   - end_on_middle_crossing: boolean to whether to end the simulation on a crossing of the
       middle milestone
  Output:
   - force: the SEEKR force object
   - data_file_name: the data file that will be monitored for crossing events.
   '''
  force = SeekrForce() # create the SEEKR force object
  neighbor1 = seekrcalc.milestones[milestone.neighbors[0]] # find the neighbor milestones
  neighbor2 = seekrcalc.milestones[milestone.neighbors[1]]
  radius1 = neighbor1.radius / 10.0 # extract neighbor milestone radii
  radius2 = milestone.radius / 10.0 # convert to nm. TODO: better way to deal with these units?
  radius3 = neighbor2.radius / 10.0
  data_file_name = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'fwd_rev', 'transition.dat') # define the file to write transition information
  # Define all settings and parameters for the SEEKR force object
  force.addSphericalMilestone(len(milestone.atom_selection_1), len(milestone.atom_selection_2), radius1, radius2, radius3, milestone.atom_selection_1, milestone.atom_selection_2, end_on_middle_crossing, data_file_name)
  system.addForce(force) # Add the SEEKR force to the openMM system
  if verbose: print "SEEKR force added to system."
  return force, data_file_name

def get_data_file_length(data_file_name):
  '''Open a file, and read the number of lines in it
  Input:
   - data_file_name: A string of the file name to read the length of
  Output:
   - file_length: An integer of the number of lines within the file
  '''
  if not os.path.exists(data_file_name): # if the file doesn't exist, return 0
    return 0
  data_file = open(data_file_name, 'r') # open the file for reading
  file_length = len(data_file.readlines()) # get the number of lines
  data_file.close()
  return file_length

def read_data_file_transitions(data_file_name, seekrcalc, milestone):
  '''Open the data file, read the transitions, and populate a dictionary of the 
  milestone transitions.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to run the simulation for
   - data_file_name: A string of the file name to read the transitions from
  Output:
   - transition_dict: Dictionary that contains transition information. Format:
       {'src_dest1':count1, 'src_dest2':count2}.
   - avg_incubation_time: Float of the average time spent in this milestone.
  '''
  transition_dict = {}
  incubation_list = []
  neighbor1 = seekrcalc.milestones[milestone.neighbors[0]] # find the neighbor milestones
  neighbor2 = seekrcalc.milestones[milestone.neighbors[1]]
  data_file = open(data_file_name, 'r')
  dest1 = neighbor1.index
  dest2 = neighbor2.index
  src = milestone.index
  key_string1 = '%d_%d' % (src, dest1)
  key_string2 = '%d_%d' % (src, dest2)
  transition_dict[key_string1] = 0
  transition_dict[key_string2] = 0
  for i, line in enumerate(data_file.readlines()):
    line = line.split()
    trans = line[0]
    time = float(line[1])
    if trans == '1': 
      transition_dict[key_string1] += 1
    elif trans == '3':
      transition_dict[key_string2] += 1
    else:
      raise AssertionError, "An unexpected value was found in the transition data file: "+str(trans)
    incubation_time_list.append(time)
  incubation_time_summation = reduce(sum, incubation_time_list)
  avg_incubation_time = incubation_time_summation / len(incubation_time_list)
  data_file.close()
  return transition_dict, avg_incubation_time

def read_data_file_successes(data_file_name):
  '''Open the data file, read the transitions, and populate a list of successful
  transitions.
  Input:
   - data_file_name: A string of the file name to read the transitions from
  Output:
   - success_list: Dictionary that contains transition information.
  '''
  success_list = []
  data_file = open(data_file_name, 'r')
  for i, line in enumerate(data_file.readlines()):
    line = line.split()
    if line[0] != '2':
      success_list.append(i)
  data_file.close()
  return success_list

def launch_fwd_rev_stage(seekrcalc, milestone, traj_base, end_on_middle_crossing, input_coords, input_vels=None, box_vectors=None):
  '''launch a reversal stage SEEKR calculation.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to run the simulation for
   - input_coords: the set of input coordinates for the trajectory
   - box_vectors: A 3x3 matrix representing the triclinic box vectors
  Output:
   - starting_positions: the set of starting positions for all simulations
   - starting_velocities: the set of starting velocities for simulations
   - data_file_name: string of the file where the transitions are stored
   - indices_list: a list of tuples that have the position and velocity indices
  '''
  # Open Amber input files
  prmtop_filename = milestone.openmm.prmtop_filename
  inpcrd_filename = milestone.openmm.inpcrd_filename
  if verbose: print "opening files:", prmtop_filename, inpcrd_filename
  prmtop = AmberPrmtopFile(prmtop_filename)
  inpcrd = AmberInpcrdFile(inpcrd_filename)
  
  # create OpenMM NVE system
  system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds) # This is fine because h-bonds are always constrained in water!
  integrator = VerletIntegrator(0.002*picoseconds)
  platform = Platform.getPlatformByName('CUDA')
  properties = seekrcalc.openmm.properties
  
  # create and prepare the SEEKR milestones
  myforce, data_file_name = create_spherical_seekr_force(seekrcalc, milestone, system, end_on_middle_crossing)
  
  os.system('rm %s' % data_file_name) # delete the existing transition data file if it exists
  
  simulation = Simulation(prmtop.topology, system, integrator, platform, properties) # create the system object
  
  starttime = time.time()
  if verbose: print "running a total of", len(input_coords) * seekrcalc.fwd_rev_stage.launches_per_config, "simulations in this stage"
  starting_velocities = [] # keep track of a list of all starting velocities
  starting_positions = [] # also keep track of starting positions
  indices_list = []
  for i in range(len(input_coords)): # for all coordinates in the list
    for j in range(seekrcalc.fwd_rev_stage.launches_per_config): # for however many times the 
      simulation.context.setPositions(input_coords[i])
      indices_list.append((i,j))
      simulation.context.setTime(0.0)
      if box_vectors:
        simulation.context.setPeriodicBoxVectors(*box_vectors)
      elif inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
        
      if input_vels == None: # if no velocities are provided, then assign by Maxwell Boltzmann
        simulation.context.setVelocitiesToTemperature(seekrcalc.master_temperature*kelvin)
      else: # assign provided velocities
        assert len(input_vels) == len(input_coords), "The lengths of input_coords and input_vels must be equal"
        simulation.context.setVelocities(input_vels[i])
      state = simulation.context.getState(getVelocities = True, getPositions = True)
      velocities = state.getVelocities()
      positions = state.getPositions()
      starting_velocities.append(deepcopy(velocities)) # deepcopy to prevent byref artifacts
      starting_positions.append(deepcopy(positions))
      traj_name = traj_base+"%d_%d.dcd" % (i, j) #"fwd_rev%d_%d.dcd" % (i, j)
      fwd_rev_traj = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'fwd_rev', traj_name)
      simulation.reporters = [StateDataReporter(stdout, seekrcalc.fwd_rev_stage.energy_freq, step=True, potentialEnergy=True, temperature=True, volume=True)]
      simulation.reporters.append(DCDReporter(fwd_rev_traj, seekrcalc.fwd_rev_stage.traj_freq))
      
      data_file_length = get_data_file_length(data_file_name)
      counter = 0
      while get_data_file_length(data_file_name) == data_file_length:
        data_file_length = get_data_file_length(data_file_name)
        simulation.step(seekrcalc.fwd_rev_stage.steps)
        counter += 1
        if counter > MAX_REVERSE_ITER: 
          print "maximum iterations exceeded."
          break
  print "Time elapsed:", time.time() - starttime
  return starting_positions, starting_velocities, data_file_name, indices_list

def process_reversal_data(reversal_coordinates, reversal_velocities, data_file_name):
  '''Reads the transition data file, then sorts the set of starting positions
  and velocities based on which ones 'succeeded'
  Input:
   - reversal coordinates: a list of all coordinates where reversals were
       started from.
   - reversal velocities: a list of all velocities where reversals were
       started from.
   - data_file_name: the file where the transitions were written
  Output:
   - success_coordinates: the coordinates, in order, where the reversals
       succeeded.
   - success_velocities: the velocities, in order, where the reversals succeeded
  '''
  success_coordinates = []
  success_velocities = []
  success_indices = read_data_file_successes(data_file_name)
  for index in success_indices:
    vel = reversal_velocities[index]
    pos = reversal_coordinates[index]
    success_coordinates.append(pos)
    success_velocities.append(vel)
  return success_coordinates, success_velocities

def pickle_coords_vels(seekrcalc, milestone, reversal_coords, reversal_vels, success_coords, success_vels):
  '''Save the reversal starting positions and velocities in a pickle for easy
  retrieval.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
   - milestone: the Milestone() object to run the simulation for
   - reversal_coords: a list of all coordinates where reversals were
       started from.
   - reversal_vels: a list of all velocities where reversals were started from
   - success_coords: a list of all coordinates where reversals succeeded
   - success_vels: a list of all coordinates where reversals succeeded
  Output:
   - reversal_coords_pickle: reversal_coords pickled
   - reversal_vels_pickle: reversal_vels pickled
   - success_coords_pickle: success_coords pickled
   - success_vels_pickle: success_vels pickled
  '''
  reversal_coords_pickle = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'fwd_rev', 'reversal_coords.pickle')
  reversal_vels_pickle = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'fwd_rev', 'reversal_vels.pickle')
  success_coords_pickle = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'fwd_rev', 'success_coords.pickle')
  success_vels_pickle = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'fwd_rev', 'success_vels.pickle')
  
  reversal_coords_pickle_file=open(reversal_coords_pickle, 'wb')
  pickle.dump(reversal_coords, reversal_coords_pickle_file, protocol=-1) # protocol=-1 means to use the fastest available protocol method
  reversal_coords_pickle_file.close()
  seekrcalc.fwd_rev_stage.reversal_coords_pickle = reversal_coords_pickle
  
  reversal_vels_pickle_file=open(reversal_vels_pickle, 'wb')
  pickle.dump(reversal_vels, reversal_vels_pickle_file, protocol=-1)
  reversal_vels_pickle_file.close()
  seekrcalc.fwd_rev_stage.reversal_vels_pickle = reversal_vels_pickle
  
  success_coords_pickle_file=open(success_coords_pickle, 'wb')
  pickle.dump(success_coords, success_coords_pickle_file, protocol=-1)
  success_coords_pickle_file.close()
  seekrcalc.fwd_rev_stage.success_coords_pickle = success_coords_pickle
  
  success_vels_pickle_file=open(success_vels_pickle, 'wb')
  pickle.dump(success_vels, success_vels_pickle_file, protocol=-1)
  success_vels_pickle_file.close()
  seekrcalc.fwd_rev_stage.success_vels_pickle = success_vels_pickle
  
  return reversal_coords_pickle, reversal_vels_pickle, success_coords_pickle, success_vels_pickle
  