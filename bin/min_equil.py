"""
This module runs minimizations and temperature equilibrations using 
OpenMM.

Created on May 10, 2018

@author: lvotapka
"""

import os
from sys import stdout

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *


verbose = True



def run_min_equil(seekrcalc):
  """Run the minimizations and temperature equilibrations of the system.
  Input:
   - seekrcalc: The SeekrCalculation object that contains all the 
       settings for the SEEKR calculation.
  Output:
   - None
  """
  for milestone in seekrcalc.milestones:
    if milestone.md and milestone.openmm.system:
      # minimize energy
      simulation = milestone.openmm.simulation
      
      if verbose: 
        print "Running energy minimization on milestone:", milestone.index
      simulation.minimizeEnergy()
      
      reporter_output = os.path.join(seekrcalc.project.rootdir, 
                                     milestone.directory, 
                                     "md", "temp_equil", "temp_equil.dcd")
      simulation.reporters.append(DCDReporter(reporter_output, 100))
      
      if verbose: print "Running temperature equilibration."
      for i, temperature in enumerate(
          seekrcalc.min_equil.temp_equil_temperatures):
        simulation.context.setVelocitiesToTemperature(temperature*kelvin)
        simulation.integrator.setTemperature(temperature*kelvin)
        simulation.step(seekrcalc.min_equil.temp_equil_steps)
        if verbose: 
          print "Ran temperature equilibration:", \
              simulation.integrator.getTemperature()
        
  return
