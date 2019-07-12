'''
Created on July 3, 2019

@author: lvotapka

Creates an OpenSEEKR calculation for the barnase-barstar system
'''

import seekr
from seekr import amber, bd
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
me.master_temperature = 298.15 #*kelvin # temperature of (most) all calculations

# project information
me.project.name = 'barnasebarstar'
me.project.rootdir = '/home/lvotapka/barnasebarstar'
me.project.empty_rootdir = remove_old_filetree
me.md = True
me.bd = False

# OpenMM information
me.openmm.platform = Platform.getPlatformByName('CUDA')
me.openmm.properties = {'CudaDeviceIndex':'0', 'CudaPrecision':'mixed'}

# Selection information
# in wet structure, SER 36 and ARG 57
barnase_indices = [833, 549] # ARG 59, SER 36
barstar_indices = [568, 687] # ASP 35, TRP 44 [2268, 2387]
me.selections.site_com_indices = barnase_indices

# Building information
me.building.ff = 'amber'
me.building.lig_dry_pqr_filename = '/home/lvotapka/barnasebarstar_files/barstar_dry.pqr'
me.building.rec_wet_pdb_filename = '/home/lvotapka/barnasebarstar_files/barnase_dry.pdb'
me.building.rec_dry_pqr_filename = '/home/lvotapka/barnasebarstar_files/barnase_dry.pqr'
me.building.reject_clashes = True

# Minimization / Temperature Equilibration info
me.min_equil.min_num_steps = 15000
me.min_equil.min_reporter_freq = 500 #[PDBReporter('dummy', 500)] # SEEKR will automatically change the filename
me.min_equil.temp_equil_integrator = LangevinIntegrator(me.master_temperature*kelvin, 5/picosecond, 0.002*picoseconds)
me.min_equil.temp_equil_reporters = [PDBReporter('dummy', 100)] # SEEKR will automatically change the filename
me.min_equil.temp_equil_steps = 10000 # number of simulation steps per temperature
me.min_equil.temp_equil_temperatures = [ 298.15] # progression of the temperature equilibration

# BrownDye information
me.browndye.browndye_bin_dir = '/home/lvotapka/Downloads/browndye/bin'
me.browndye.num_threads = 10
me.browndye.lig_dry_pqr_filename = '/home/lvotapka/barnasebarstar_files/barstar_dry.pqr'
me.browndye.rec_dry_pqr_filename = '/home/lvotapka/barnasebarstar_files/barnase_dry.pqr'
me.browndye.prods_per_anchor = 1000000
me.browndye.apbs.executable = '/home/lvotapka/Downloads/APBS-1.5-linux64/bin/apbs'
me.browndye.fhpd_numtraj = 1000

ion1 = seekr.APBS_ion('Cl-', 0.05, -1.0, 1.67)
ion2 = seekr.APBS_ion('Na+', 0.05, 1.0, 1.16) # define ions in the system


me.browndye.apbs.ions = [ion1, ion2]
me.browndye.apbs.linear_pbe = False
me.browndye.apbs.inputgen.executable = '/home/lvotapka/Downloads/APBS-1.5-linux64/share/apbs/tools/manip/inputgen.py'

# Generate Milestones
#radius_list = np.arange(2.0, 14.1, 2.0)
radius_list = [8.0, 12.0, 16.0, 20.0, 24.0, 28.0]
vector = np.array([6.926, -5.359, -5.908])
barnase_points = [np.array([24.872, 37.120, -1.407]), np.array([39.597, 42.427, 15.375])]
barstar_points = [np.array([24.626, 36.903, 2.106]), np.array([35.476, 39.775, 12.871])]
milestones = seekr.generate_RMSD_milestones(me, barnase_indices, barstar_indices, barnase_points, barstar_points, radius_list, 0, vector, absolute=False)
print "The following milestones were created:"
seekr.print_RMSD_milestone_info(milestones)
me.milestones = milestones

# Generate Filetree and Building files
seekr.generate_filetree(me)
holo_config_wet, insert_index, last_ligand_index = seekr.generate_configs(me)
ligand_heavy_indices = seekr.find_heavy_atoms(holo_config_wet, range(insert_index, insert_index+last_ligand_index+1))
print "Insert_index:", insert_index
print "Ligand heavy indices:", ligand_heavy_indices  # ATTN: changed

amber_settings = amber.AmberSettings()
amber_settings.leap_program = 'tleap'

amber_settings.leap_template = '''
source leaprc.protein.ff14SB
source leaprc.gaff
source leaprc.water.tip4pew
set default PBRadii mbondi2

WAT = T4E
HOH = T4E

loadAmberParams frcmod.ionsjc_tip4pew
loadAmberParams frcmod.tip4pew

barnase_barstar = loadpdb $HOLO_WET_PDB

solvateoct barnase_barstar TIP4PEWBOX 8

addIons2 barnase_barstar Na+ 0
addIons2 barnase_barstar Cl- 0

addIons2 barnase_barstar Na+ 13
addIons2 barnase_barstar Cl- 13

saveamberparm barnase_barstar $PRMTOP $INPCRD
savepdb barnase_barstar  $LEAP_OUTPUT_PDB
check barnase_barstar
charge barnase_barstar
quit

'''

me.min_equil.constrained += range(insert_index)
me.min_equil.constrained += ligand_heavy_indices
    
for milestone in me.milestones:
  milestone.atom_selection_1 = me.selections.site_com_indices
  milestone.atom_selection_2 = ligand_heavy_indices
  milestone.atom_indices2 = map(lambda x: x+insert_index, milestone.atom_indices2)
  if milestone.md:
    amber.amber_building(me, milestone, amber_settings)
    
    if not milestone.openmm.prmtop_filename: continue
    # modify the file to have the correct solvent octahedron box
    '''
    parm = pmd.load_file(milestone.openmm.prmtop_filename, xyz=milestone.openmm.inpcrd_filename)
    
    #TODO: straighten out these units
    
    parm.box = np.array([68.43194878035462, 68.43194878035462, 68.43194878035462,109.471219, 109.471219,109.471219])
    box_vector = Quantity([[68.43194878035462, 0.0, 0.0], [-22.81064775292344, 64.51826069392389, 0.0], [-22.81064775292344, -32.259126442612187, 55.87445502310258]], unit=angstrom)
    parm.box_vectors = box_vector
    print "parm.box:", parm.box
    print "parm.box_vectors:", parm.box_vectors
    print "saving prmtop for milestone:", milestone.index
    parm.save(milestone.openmm.prmtop_filename, overwrite = True)
    print "saving inpcrd for milestone:", milestone.index
    parm.save(milestone.openmm.inpcrd_filename, overwrite = True)
    '''
    amber.create_simulation(me, milestone)
    #milestone.openmm.simulation.context.setPeriodicBoxVectors([6.843194878035462, 0.0, 0.0], [-2.281064775292344, 6.451826069392389, 0.0], [-2.281064775292344, -3.2259126442612187, 5.587445502310258])
    
    if not me.openmm.simulation: # create a sample of a simulation file for future use
      me.openmm.simulation = milestone.openmm.simulation
    
seekr.run_min_equil(me)

# save equilibration output
for milestone in me.milestones:
  if milestone.md:
    filename = amber.save_restart(me, milestone)
    milestone.openmm.umbrella_pdb_filename = filename

print "Saving all system settings for Umbrella stage and later stages."
me.save()

if bd:
  print "Preparing Brownian dynamics stages..."
  bd.build_bd(me)


print "Ligand heavy indices:", ligand_heavy_indices
print "DON'T TRUST THE FOLLOWING NUMBERS"
for milestone in me.milestones:
  print "milestone[%d].atom_indices2:" % milestone.index, milestone.atom_indices2

