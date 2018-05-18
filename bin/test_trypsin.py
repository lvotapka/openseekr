'''
Created on May 9, 2018

@author: lvotapka

This is a module to test the SEEKR package and to aid in its development and
introduction of functionality
'''

import seekr
from seekr import amber
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import numpy as np
import parmed as pmd
import sys

remove_old_filetree = False
if 'remove' in sys.argv[1:]:
  remove_old_filetree = True

# Define settings object for all simulations

me = seekr.SeekrCalculation() # create a new SEEKR calculation object
me.master_temperature = 300. # temperature of (most) all calculations

# project information
me.project.name = 'test_tryp'
me.project.rootdir = '/home/lvotapka/tryp_test'
me.project.k_off = True # make sure to do a k-off calculation
me.project.empty_rootdir = remove_old_filetree

# OpenMM information
me.openmm.platform = Platform.getPlatformByName('CUDA')
me.openmm.properties = {'CudaDeviceIndex':'0', 'CudaPrecision':'single'}

# Selection information
me.selections.lig_com_indices = range(3222, 3240) # range of ligand atom indices
me.selections.site_com_indices = [2467, 2479, 2490, 2536, 2746, 2770, 2788]

# Building information
me.building.ff = 'amber'
me.building.lig_dry_pqr_filename = '/home/lvotapka/tryp_files/benzamidine_moved.pqr'
me.building.rec_wet_pdb_filename = '/home/lvotapka/openseekr/openmm_examples/leaptest/apo_postleap2.pdb'
me.building.rec_dry_pqr_filename = '/home/lvotapka/tryp_files/tryp_dry.pqr'
me.building.reject_clashes = True

# Minimization / Temperature Equilibration info
me.min_equil.min_constrained = [] #[2467, 2479, 2490, 2536, 2746, 2770, 2788] + range(3222, 3240)
me.min_equil.min_num_steps = 5000
me.min_equil.min_reporter_freq = 500 #[PDBReporter('dummy', 500)] # SEEKR will automatically change the filename
me.min_equil.temp_equil_constrained = [] # [2467, 2479, 2490, 2536, 2746, 2770, 2788] + range(3222, 3240)
me.min_equil.temp_equil_integrator = LangevinIntegrator(me.master_temperature*kelvin, 1/picosecond, 0.002*picoseconds)
me.min_equil.temp_equil_reporters = [PDBReporter('dummy', 100)] # SEEKR will automatically change the filename
me.min_equil.temp_equil_steps = 1000 # number of simulation steps per temperature
me.min_equil.temp_equil_temperatures = [300., 310., 320., 330., 340., 350., 340., 330., 320., 310., 300.] # progression of the temperature equilibration

# Umbrella info
me.umbrella_stage.integrator = LangevinIntegrator(me.master_temperature*kelvin, 1/picosecond, 0.002*picoseconds)
me.umbrella_stage.reporters = [DCDReporter('dummy', 1000)] # SEEKR will automatically change the filename
me.umbrella_stage.steps = 1000000 # number of steps to run in umbrella stage

# Fwd-rev info

me.fwd_rev_stage.integrator = VerletIntegrator(2/femtoseconds) # fwd-rev stage must be an NVE integrator
me.fwd_rev_stage.reporters = [PDBReporter('dummy', 1000)] # SEEKR will automatically change the filename
me.fwd_rev_stage.launches_per_config = 10 # number of times to reinitialize velocity to start in the fwd-rev stage

# BrownDye information
me.browndye.browndye_bin_dir = '/home/lvotapka/Downloads/browndye/bin'
me.browndye.num_threads = 10
me.browndye.lig_dry_pqr_filename = '/home/lvotapka/tryp_files/benzamidine_moved.pqr'
me.browndye.rec_dry_pqr_filename = '/home/lvotapka/tryp_files/tryp_dry.pqr'
me.browndye.prods_per_anchor = 1000000
me.browndye.apbs.executable = '/home/lvotapka/Downloads/APBS-1.5-linux64/apbs'

ion1 = seekr.APBS_ion('Cl-', 0.10, -1.0, 1.67)
ion2 = seekr.APBS_ion('Ca2+', 0.02, 2.0, 2.1) # define ions in the system
ion3 = seekr.APBS_ion('tris', 0.06, 1.0, 4.0)

me.browndye.apbs.ions = [ion1, ion2, ion3]
me.browndye.apbs.linear_pbe = False
me.browndye.apbs.inputgen.executable = '/home/lvotapka/Downloads/APBS-1.5-linux64/share/apbs/tools/manip'

# Generate Milestones
rec_site_atom_indices = [2467, 2479, 2490, 2536, 2746, 2770, 2788]
origin = np.array([-2.777, 11.117, 0.397])
radius_list = np.arange(2.0, 14.1, 2.0)
vectors = [np.array([1.0, 1.0, 1.0]), np.array([-3.0, 0.0, 4.0,])]
milestones = seekr.generate_spherical_milestones(rec_site_atom_indices, origin, radius_list, 0, vectors, k_off=False, absolute=False)
print "The following milestones were created:"
seekr.print_spherical_milestone_info(milestones)
me.milestones = milestones

# Generate Filetree and Building files
seekr.generate_filetree(me)
seekr.generate_configs(me)

amber_settings = amber.AmberSettings()
amber_settings.leap_program = 'tleap'

amber_settings.leap_template = '''
source leaprc.protein.ff14SB
source leaprc.ff14SB.redq
source leaprc.gaff
set default FlexibleWater on
set default PBRadii mbondi2
loadoff /home/lvotapka/torq/lvotapka/Documents/trypsin_files/Ca2.lib
loadoff /home/lvotapka/torq/lvotapka/Documents/trypsin_files/benzamidine.lib
loadamberparams /home/lvotapka/torq/lvotapka/Documents/trypsin_files/benzamidine.frcmod
WAT = T4E
HOH = T4E
loadAmberParams frcmod.ionsjc_tip4pew
loadAmberParams frcmod.tip4pew


holo = loadpdb $HOLO_WET_PDB
bond holo.7.SG holo.137.SG
bond holo.25.SG holo.41.SG
bond holo.109.SG holo.210.SG
bond holo.116.SG holo.183.SG
bond holo.148.SG holo.162.SG
bond holo.173.SG holo.197.SG

addIons2 holo Cl- 0


charge holo
check holo

saveamberparm holo $PRMTOP $INPCRD

savepdb holo  $LEAP_OUTPUT_PDB

quit
'''

me.min_equil.constrained += range(3238)
    
for milestone in me.milestones:
  milestone.atom_selection_1 = me.selections.site_com_indices
  milestone.atom_selection_2 = me.selections.lig_com_indices
  if milestone.md:
    amber.amber_building(me, milestone, amber_settings)
    
    if not milestone.openmm.prmtop_filename: continue
    # modify the file to have the correct solvent octahedron box
    #print "modifying prmtop/inpcrd pair to have the correct solvent box"
    parm = pmd.load_file(milestone.openmm.prmtop_filename, xyz=milestone.openmm.inpcrd_filename)
    
    parm.box = np.array([64.9127105, 64.9127105, 64.9127105,109.471219, 109.471219,109.471219])
    box_vector = Quantity([[64.912710500000003, 0.0, 0.0], [-21.637568420791037, 61.200290990259163, 0.0], [-21.637568420791037, -30.600141791568205, 53.00100885481632]], unit=angstrom)
    parm.box_vectors = box_vector
    print "parm.box:", parm.box
    print "parm.box_vectors:", parm.box_vectors
    print "saving prmtop for milestone:", milestone.index
    parm.save(milestone.openmm.prmtop_filename, overwrite = True)
    print "saving inpcrd for milestone:", milestone.index
    parm.save(milestone.openmm.inpcrd_filename, overwrite = True)
    
    amber.create_simulation(me, milestone)
    milestone.openmm.simulation.context.setPeriodicBoxVectors([6.4912710500000003, 0.0, 0.0], [-2.1637568420791037, 6.1200290990259163, 0.0], [-2.1637568420791037, -3.0600141791568205, 5.300100885481632])
    
    if not me.openmm.simulation: # create a sample of a simulation file for future use
      me.openmm.simulation = milestone.openmm.simulation


    
# make ligand restrained
me.min_equil.constrained += me.selections.lig_com_indices
    
seekr.run_min_equil(me)

# save equilibration output
me.openmm
for milestone in me.milestones:
  if milestone.md:
    filename = amber.save_restart(me, milestone)
    milestone.openmm.umbrella_pdb_filename = filename

print "Saving all system settings for Umbrella stage and later stages."
me.save()