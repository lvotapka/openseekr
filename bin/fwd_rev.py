"""
Functions and objects intended to run the reversal and forward stages of 
a SEEKR calculation.

Created on June 15, 2018

@author: lvotapka
"""
import sys
import os
import time
from sys import stdout
from copy import deepcopy
import cPickle as pickle


import parmed
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from seekr import amber
from seekrplugin import SeekrForce
from adv_template import Adv_template



verbose = True

MAX_REVERSE_ITER = 100000 # 200 ns
MAX_FORWARD_ITER = 100000 # 200 ns

def make_box_info(box_vectors):
  """Converts a 3x3 matrix of box vectors to x,y,z,alpha,beta,gamma 
  format for CPPTRAJ input.
  Input:
   - box_vectors: A 3x3 matrix representing the triclinic box vectors
  Output:
   - box_info_string: a string formatted for CPPTRAJ 'box' command
  """

  # extract the number without units
  box_length = box_vectors[0][0].value_in_unit(angstroms) 
  box_info_string = "x %f y %f z %f alpha 109.4712190 beta 109.4712190  \
      gamma 109.4712190" % (box_length, box_length, box_length)
  return box_info_string

def autoimage_traj(parm_name, trajin_name, trajout_name, box_info, 
                   cpptraj_script_location, cpptraj_exe="cpptraj", 
                   writing_frames=()):
  """Runs the CPPTRAJ autoimage command for a triclinic box simulation.
  Input:
   - parm_name: a string representing the filename of a .prmtop or 
       .parm7AM BER parameter/topology file
   - trajin_name: a non-imaged trajectory to load for imaging
   - box_info: a string representing the triclinic box in x,y,z,alpha,
        beta,gamma format
   - cpptraj_script_location: string for the location to write the 
        cpptraj script
   - cpptraj_exe: an optional string representing the OS command to 
        run CPPTRAJ
  Output:
   - None
   """
  cpptraj_template = """parm $PARMFILE
box $BOX_INFO
trajin $TRAJIN
autoimage
trajout $TRAJOUT $FRAME_STR
go
quit
"""
  assert len(writing_frames) < 4, "When writing the autoimaged trajectory, the\
      format of the writing_frames variable must be: (start, stop, offset)"
  if len(writing_frames) == 0:
    # then write all frames
    frame_str = ""
  elif len(writing_frames) == 1:
    # only include the start
    frame_str = "start %d" % writing_frames[0]
  elif len(writing_frames) == 2:
    frame_str = "start %d stop %d" % (writing_frames[0], writing_frames[1])
  else: # the length is 3
    frame_str = "start %d stop %d offset %d" % (writing_frames[0], 
                                                writing_frames[1], 
                                                writing_frames[2])
  # define template dictionary
  cpptraj_dict = {"PARMFILE":parm_name, "TRAJIN":trajin_name, 
                  "TRAJOUT":trajout_name, "BOX_INFO":box_info, 
                  "FRAME_STR":frame_str}
  # fill in the values into the template from the dictionary
  cpptraj_script = Adv_template(cpptraj_template, cpptraj_dict)
  # open the script for writing
  extract_file = open(cpptraj_script_location, "w")
  # write a cpptraj script
  extract_file.write(cpptraj_script.get_output())
  extract_file.close()
  # run cpptraj
  os.system("%s < %s" % (cpptraj_exe, cpptraj_script_location))
  return

def create_spherical_seekr_force(seekrcalc, milestone, system, 
                                 end_on_middle_crossing, 
                                 transition_filename="transition.dat"):
  """create the SEEKR 'force' even though it's more of a monitor than a 
  force.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the 
        settings for the SEEKR calculation.
   - milestone: the Milestone() object to run the simulation for
   - system: the OpenMM system to add the SEEKR 'force' to
   - end_on_middle_crossing: boolean to whether to end the simulation 
        on a crossing of the middle milestone
  Output:
   - force: the SEEKR force object
   - data_file_name: the data file that will be monitored for 
        crossing events.
   """
  # create the SEEKR force object
  force = SeekrForce()
  if len(milestone.neighbors) == 2:
    # find the neighbor milestones
    neighbor1 = seekrcalc.milestones[milestone.neighbors[0]]
    neighbor2 = seekrcalc.milestones[milestone.neighbors[1]]
    # extract neighbor milestone radii
    radius1 = neighbor1.radius / 10.0
    # convert to nm. TODO: better way to deal with these units?
    radius2 = milestone.radius / 10.0
    radius3 = neighbor2.radius / 10.0
  # this is an endpoint milestone HACKY...
  elif len(milestone.neighbors) == 1:
    # find the neighbor milestones
    neighbor1 = seekrcalc.milestones[milestone.neighbors[0]]
    neighbor2 = seekrcalc.milestones[milestone.neighbors[0]]
    radius1 = 0.0
    # extract neighbor milestone radii
    # convert to nm. TODO: better way to deal with these units?
    radius2 = milestone.radius / 10.0
    radius3 = neighbor2.radius / 10.0
  else:
    raise Exception, "Only one or two milestone neighbors allowed at present. \
        Number of neighbors: %d" % milestone.neighbors
  # define the file to write transition information
  data_file_name = os.path.join(seekrcalc.project.rootdir, 
                                milestone.directory,
                                "md", "fwd_rev", transition_filename) 
  # Define all settings and parameters for the SEEKR force object
  force.addSphericalMilestone(len(milestone.atom_selection_1), 
                              len(milestone.atom_selection_2), 
                              radius1, radius2, radius3, 
                              milestone.atom_selection_1, 
                              milestone.atom_selection_2, 
                              end_on_middle_crossing, data_file_name)
  # Add the SEEKR force to the openMM system
  system.addForce(force)
  if verbose: print "SEEKR force added to system."
  return force, data_file_name

def get_data_file_length(data_file_name):
  """Open a file, and read the number of lines in it
  Input:
   - data_file_name: A string of the file name to read the length of
  Output:
   - file_length: An integer of the number of lines within the file
  """
  # if the file doesn't exist, return 0
  if not os.path.exists(data_file_name):
    return 0
  # open the file for reading
  data_file = open(data_file_name, "r")
  # get the number of lines
  file_length = len(data_file.readlines())
  data_file.close()
  return file_length

def read_data_file_transitions(data_file_name, seekrcalc, milestone):
  """Open the data file, read the transitions, and populate a dictionary 
  of the milestone transitions.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the 
       settings for the SEEKR calculation.
   - milestone: the Milestone() object to run the simulation for
   - data_file_name: A string of the file name to read the 
       transitions from
  Output:
   - transition_dict: Dictionary that contains transition information. 
       Format:{'src_dest1':count1, 'src_dest2':count2}.
   - avg_incubation_time: Float of the average time spent in this 
       milestone.
  """
  transition_dict = {}
  incubation_time_list = []
  num_failed = 0
  data_file = open(data_file_name, "r")
  if len(milestone.neighbors) == 2:
    # find the neighbor milestones
    neighbor1 = seekrcalc.milestones[milestone.neighbors[0]]
    neighbor2 = seekrcalc.milestones[milestone.neighbors[1]]
    dest1 = neighbor1.index
    dest2 = neighbor2.index
  else:
    neighbor2 = seekrcalc.milestones[milestone.neighbors[0]]
    dest2 = neighbor2.index
    dest1 = -1
  src = milestone.index
  key_string1 = "%d_%d" % (src, dest1)
  key_string2 = "%d_%d" % (src, dest2)
  transition_dict[key_string1] = 0
  transition_dict[key_string2] = 0
  for i, line in enumerate(data_file.readlines()):
    line = line.split()
    trans = line[0]
    time = float(line[1])
    if trans == "1": 
      transition_dict[key_string1] += 1
    elif trans == "3":
      transition_dict[key_string2] += 1
    elif trans == "3*" or trans == "1*":
      num_failed += 1
    else:
      raise AssertionError, "An unexpected value was found in the transition \
          data file: "+str(trans)
    incubation_time_list.append(time)
  incubation_time_summation = sum(incubation_time_list)
  avg_incubation_time = incubation_time_summation / len(incubation_time_list)
  data_file.close()
  print "number of failed forward trajectories:", num_failed
  return transition_dict, avg_incubation_time, incubation_time_list

def read_data_file_successes(data_file_name):
  """Open the data file, read the transitions, and populate a list of 
  successful transitions.
  Input:
   - data_file_name: A string of the file name to read the transitions 
       from
  Output:
   - success_list: Dictionary that contains transition information.
  """
  success_list = []
  data_file = open(data_file_name, "r")
  for i, line in enumerate(data_file.readlines()):
    line = line.split()
    if line[0] != "2":
      success_list.append(i)
  data_file.close()
  return success_list

def read_reversal_data_file_last(data_file_name):
  """Open the data file, read the final transition, and return whether 
  the latest reversal was successful.
  Input:
   - data_file_name: A string of the file name to read the transitions 
       from
  Output:
   - success: A boolean about whether the latest trajectory was 
       successful
  """
  data_file = open(data_file_name, "r")
  line = data_file.readlines()[-1]
  data_file.close()
  if line[0] != "2":
    return True
  else:
    return False
  
def sort_pickle_key(pickle_filename):
  """A function used by the 'key' argument of the sorted() function to 
  sort the forward dcd files."""
  # TODO: this is hacky..
  basenum = os.path.basename(pickle_filename).split(".")[0][14:] .
  sorting_number = int(basenum)
  return sorting_number

def sort_forward_dcd_key(forward_dcd):
  """A function used by the 'key' argument of the sorted() function 
  to sort the forward dcd files.
  """
  # TODO: this is hacky...
  basename = os.path.basename(forward_dcd).split(".")[0][7:]
  index_list = basename.split("_")
  if len(index_list) == 3:
    sorting_list = [int(index_list[2]), int(index_list[0]), int(index_list[1])]
  elif len(index_list) == 2: # backwards compatibility
    sorting_list = [int(index_list[0]), int(index_list[1])]
  return sorting_list

def launch_fwd_rev_stage(seekrcalc, milestone, traj_base, 
                         end_on_middle_crossing, dcd_iterator, 
                         dcd_iterator_chunk=9e9, input_vels=None, 
                         box_vectors=None, 
                         transition_filename="transition.dat", suffix="", 
                         save_fwd_rev=False, save_last_frame=True):
  # TODO: update the docstring to current inputs
  """launch a reversal stage SEEKR calculation.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the 
       settings for the SEEKR calculation.
   - milestone: the Milestone() object to run the simulation for
   - input_coords: the set of input coordinates for the trajectory
   - box_vectors: A 3x3 matrix representing the triclinic box vectors
  Output:
   - starting_positions: the set of starting positions for all 
       simulations
   - starting_velocities: the set of starting velocities for simulations
   - data_file_name: string of the file where the transitions are stored
   - indices_list: a list of tuples that have the position and 
       velocity indices
  """
  # Open Amber input files
  # then it's a forward stage, but there shouldn't be any 
  # velocity randomizations
  if input_vels != None and seekrcalc.fwd_rev_stage.launches_per_config != 1:
    raise Exception, "The variable 'fwd_rev_stage.launches_per_config' should \
        be set to 1 for the forward stage. (or velocities should not be \
        provided to this function)."
  
  prmtop_filename = milestone.openmm.prmtop_filename
  inpcrd_filename = milestone.openmm.inpcrd_filename
  if verbose: print "opening files:", prmtop_filename, inpcrd_filename
  prmtop = AmberPrmtopFile(prmtop_filename)
  inpcrd = AmberInpcrdFile(inpcrd_filename)
  
  # create OpenMM NVE system
  # This is fine because h-bonds are always constrained in water!
  system = prmtop.createSystem(nonbondedMethod=PME, 
                               nonbondedCutoff=1*nanometer, 
                               constraints=HBonds)
  integrator = VerletIntegrator(0.002*picoseconds)
  platform = Platform.getPlatformByName("CUDA")
  properties = seekrcalc.openmm.properties
  
  # create and prepare the SEEKR milestones
  myforce, data_file_name = create_spherical_seekr_force(
      seekrcalc, milestone, system, 
      end_on_middle_crossing, transition_filename)
  
  if save_fwd_rev == False:
    print "Deleting transition file:", data_file_name
    # delete the existing transition data file if it exists
    os.system("rm %s" % data_file_name)
  
  # create the system object
  simulation = Simulation(prmtop.topology, system, 
                          integrator, platform, properties)
  
  starttime = time.time()
  success_velocities = []
  success_positions = []
  indices_list = []
  i = 0
  num_errors = 0
  complete = True
  
  for dcd_frame in dcd_iterator:
    # for however many times the
    for j in range(seekrcalc.fwd_rev_stage.launches_per_config): 
      had_error = False
      if hasattr(dcd_frame, "xyz"):
        simulation.context.setPositions(dcd_frame.xyz[0])
      else:
        simulation.context.setPositions(dcd_frame)
      indices_list.append((i,j))
      simulation.context.setTime(0.0)
      if box_vectors:
        simulation.context.setPeriodicBoxVectors(*box_vectors)
      elif inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
      
      # if no velocities are provided, then assign by Maxwell Boltzmann
      if input_vels == None:
        print "assigning random velocities"
        simulation.context.setVelocitiesToTemperature(
            seekrcalc.master_temperature*kelvin)
      else: # assign provided velocities
        print "assigning provided velocities"
        simulation.context.setVelocities(input_vels[i])
      state = simulation.context.getState(getVelocities = True, 
                                          getPositions = True)
      velocities = state.getVelocities()
      positions = state.getPositions()
      #"fwd_rev%d_%d.dcd" % (i, j)
      traj_name = traj_base+"%d_%d%s.dcd" % (i, j, suffix)
      fwd_rev_traj = os.path.join(seekrcalc.project.rootdir, 
                                  milestone.directory, "md", 
                                  "fwd_rev", traj_name)
      simulation.reporters = [StateDataReporter(
          stdout, seekrcalc.fwd_rev_stage.energy_freq, step=True, 
          potentialEnergy=True, temperature=True, volume=True)]
      simulation.reporters.append(DCDReporter(fwd_rev_traj, 
                                  seekrcalc.fwd_rev_stage.traj_freq))
      
      data_file_length = get_data_file_length(data_file_name)
      counter = 0
      while get_data_file_length(data_file_name) == data_file_length:
        data_file_length = get_data_file_length(data_file_name)
        try:
          simulation.step(seekrcalc.fwd_rev_stage.steps)
        # if there was a NAN error
        except Exception:
          print "Error encountered. Continuing with the next frame."
          num_errors += 1
          had_error = True
          # don't want to log this as a success
          break
          
        counter += 1
        if counter > MAX_REVERSE_ITER: 
          print "maximum iterations exceeded."
          break
      if had_error == True:
        # move on to the next frame
        break
      if read_reversal_data_file_last(data_file_name):
        milestone.openmm.simulation = simulation
        success_positions.append(positions)
        success_velocities.append(velocities)
        if save_last_frame:
          #"fwd_rev%d_%d.dcd" % (i, j)
          pdb_last_frame_base_name = traj_base+"%d_%d%s.pdb" % (i, j, suffix)
          pdb_last_frame_name = os.path.join(seekrcalc.project.rootdir, 
                                             milestone.directory, 
                                             "md", "fwd_rev", 
                                             pdb_last_frame_base_name)
          if seekrcalc.building.ff == "amber":
            amber_parm = parmed.amber.AmberParm(
                milestone.openmm.prmtop_filename, 
                milestone.openmm.inpcrd_filename)
            amber_parm.positions = positions
            amber_parm.save(pdb_last_frame_name, overwrite=True)
          elif seekrcalc.building.ff == "charmm":
            raise Exception, "charmm ff not yet implemented."
          else:
            raise Exception, "ff not yet implemented: %s" \
                % seekrcalc.building.ff
        
    i += 1
    
    if len(success_positions) >= dcd_iterator_chunk:
      complete = False
      break
    #if i*seekrcalc.fwd_rev_stage.launches_per_config >= \
    #   dcd_iterator_chunk:
    #  break
      
  print "Time elapsed:", time.time() - starttime
  print "Number of errors:", num_errors
  
  return success_positions, success_velocities, \
      data_file_name, indices_list, complete

def process_reversal_data(reversal_coordinates, 
                          reversal_velocities, 
                          data_file_name):
  """Reads the transition data file, then sorts the set of starting 
  positions and velocities based on which ones 'succeeded'
  Input:
   - reversal coordinates: a list of all coordinates where reversals 
       were started from.
   - reversal velocities: a list of all velocities where reversals were
       started from.
   - data_file_name: the file where the transitions were written
  Output:
   - success_coordinates: the coordinates, in order, where the reversals
       succeeded.
   - success_velocities: the velocities, in order, where the reversals 
       succeeded
  """
  success_coordinates = []
  success_velocities = []
  success_indices = read_data_file_successes(data_file_name)
  for index in success_indices:
    vel = reversal_velocities[index]
    pos = reversal_coordinates[index]
    success_coordinates.append(pos)
    success_velocities.append(vel)
  return success_coordinates, success_velocities

def pickle_coords_vels(seekrcalc, milestone, 
                       success_coords, success_vels, index=None):
  """Save the reversal starting positions and velocities in a pickle 
  for easy retrieval.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the 
       settings for the SEEKR calculation.
   - milestone: the Milestone() object to process for
   - success_positions: a list of all coordinates where reversals 
       succeeded
   - success_velocities: a list of all coordinates where reversals 
       succeeded
  Output:
   - success_coords_pickle: success_coords pickled
   - success_vels_pickle: success_vels pickled
  """
  if index == None:
    suffix = ""
  else:
    suffix = str(index)
  success_coords_pickle = os.path.join(seekrcalc.project.rootdir, 
                                       milestone.directory, "md", "fwd_rev", 
                                       "success_coords%s.pickle" % suffix)
  success_vels_pickle = os.path.join(seekrcalc.project.rootdir, 
                                     milestone.directory, "md", "fwd_rev", 
                                     "success_vels%s.pickle" % suffix)
  
  success_coords_pickle_file=open(success_coords_pickle, "wb")
  pickle.dump(success_coords, success_coords_pickle_file, protocol=-1)
  success_coords_pickle_file.close()
  seekrcalc.fwd_rev_stage.success_coords_pickle = success_coords_pickle
  
  success_vels_pickle_file=open(success_vels_pickle, "wb")
  pickle.dump(success_vels, success_vels_pickle_file, protocol=-1)
  success_vels_pickle_file.close()
  seekrcalc.fwd_rev_stage.success_vels_pickle = success_vels_pickle
  
  return success_coords_pickle, success_vels_pickle

def pickle_transition_info(seekrcalc, milestone, 
                           transition_dict, avg_incubation_time):
  """Save the reversal starting positions and velocities in a pickle 
  for easy retrieval.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the 
       settings for the SEEKR calculation.
   - milestone: the Milestone() object to process for
   - transition_dict: a dictionary of milestone transitions, where the 
       keys are milestone indices and describe counts of transitions 
       to adjacent milestones
   - avg_incubation_time: The float representing the average time spent 
       after crossing this milestone, and before crossing another.
  Output:
   - None
  """
  transition_info_pickle = os.path.join(seekrcalc.project.rootdir, 
                                        milestone.directory, "md", 
                                        "fwd_rev", "transition_info.pickle")
  
  transition_info_pickle_file=open(transition_info_pickle,"wb")
  # NOTE: This will require two load calls to remove both 
  # objects from this pickle
  pickle.dump(transition_dict, transition_info_pickle_file, protocol=-1)
  pickle.dump(avg_incubation_time, transition_info_pickle_file, protocol=-1) 
  transition_info_pickle_file.close()
  
  return
