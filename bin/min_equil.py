'''
This module runs minimizations and temperature equilibrations using OpenMM.

Created on May 10, 2018

@author: lvotapka
'''

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import amber

import os

verbose = True



def run_min_equil(seekrcalc):
  '''Run the minimizations and temperature equilibrations of the system.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the settings for 
       the SEEKR calculation.
  Output:
   - None
  '''
  for milestone in seekrcalc.milestones:
    if milestone.md and milestone.openmm.system:
      # minimize energy
      simulation = milestone.openmm.simulation
      
      state = simulation.context.getState()
      print "box_vectors:"
      print(state.getPeriodicBoxVectors())
      
      if verbose: print "Running energy minimization on milestone:", milestone.index
      
      
      reporter_output = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'temp_equil', 'temp_equil.dcd')
      '''
      if os.path.exists(reporter_output):
        print "Temperature equilibration trajectories found. Skipping minimization and temperature equilibration."
        return
      '''
      simulation.minimizeEnergy()
      after_min_structure_name = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'min', 'after_min.pdb')
      amber.save_restart(seekrcalc, milestone, pdb_save_filename=after_min_structure_name)
      simulation.reporters.append(DCDReporter(reporter_output, 1000)) # TODO: allow users to modify this quantity
      
      if verbose: print "Running temperature equilibration."
      for i, temperature in enumerate(seekrcalc.min_equil.temp_equil_temperatures):
        simulation.context.setVelocitiesToTemperature(temperature*kelvin)
        simulation.integrator.setTemperature(temperature*kelvin)
        
        box_vectors = simulation.context.getState().getPeriodicBoxVectors()
        milestone.box_vectors = box_vectors
        print "simulation_box_vectors:", box_vectors
        
        simulation.step(seekrcalc.min_equil.temp_equil_steps)
        if verbose: print "Ran temperature equilibration:", simulation.integrator.getTemperature()
        
  return
