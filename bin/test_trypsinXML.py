'''
Created on March 24, 2020

@author: lvotapka

Test the new XML serializer
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
rebuild_amber = False
if 'remove' in sys.argv[1:]:
    remove_old_filetree = True

if 'amber' in sys.argv[1:]:
    rebuild_amber = True

# Define settings object for all simulations
me = seekr.SeekrCalculation() # create a new SEEKR calculation object
me.master_temperature = 298. #*kelvin # temperature of (most) all calculations

# project information
me.project.name = 'test_trypXML'
me.project.rootdir = '/home/lvotapka/tryp_testXML'
me.project.empty_rootdir = remove_old_filetree
me.md = True
me.bd = True

# OpenMM information
me.openmm.platform = Platform.getPlatformByName('Reference')
me.openmm.properties = {} #{'CudaDeviceIndex':'1', 'CudaPrecision':'mixed'}

# Selection information
rec_site_atom_indices = [2478, 2489, 2499, 2535, 2718, 2745, 2769, 2787, 2794, 2867, 2926] # Make sure this is selected by "index" in VMD
#2478 2489 2499 2535 2718 2745 2769 2787 2794 2867 2926
me.selections.site_com_indices = rec_site_atom_indices

# Building information
me.building.ff = 'amber'
me.building.lig_dry_pqr_filename = '/home/lvotapka/tryp_files/benzamidine_moved.pqr'
me.building.rec_wet_pdb_filename = '/home/lvotapka/torq/lvotapka/Documents/trypsin_files/apo_redo/tryp_wet_1946.pdb'
me.building.rec_dry_pqr_filename = '/home/lvotapka/torq/lvotapka/Documents/trypsin_files/apo_redo/tryp_dry_1946.pqr'
me.building.reject_clashes = True

# Minimization / Temperature Equilibration info
#me.min_equil.constrained += list(range(3220))
me.min_equil.min_num_steps = 5000
me.min_equil.min_reporter_freq = 500 #[PDBReporter('dummy', 500)] # SEEKR will automatically change the filename
me.min_equil.temp_equil_integrator = LangevinIntegrator(me.master_temperature*kelvin, 5/picosecond, 0.002*picoseconds)
me.min_equil.temp_equil_reporters = [PDBReporter('dummy', 100)] # SEEKR will automatically change the filename
me.min_equil.temp_equil_steps = 1000 # number of simulation steps per temperature
me.min_equil.temp_equil_temperatures = [300., 310., 300.] # progression of the temperature equilibration

# BrownDye information
me.browndye.browndye_bin_dir = '/home/lvotapka/Downloads/browndye/bin'
me.browndye.num_threads = 10
me.browndye.lig_dry_pqr_filename = '/home/lvotapka/tryp_files/benzamidine_moved.pqr'
me.browndye.rec_dry_pqr_filename = '/home/lvotapka/tryp_files/tryp_dry.pqr'
me.browndye.prods_per_anchor = 1000000
me.browndye.apbs.executable = '/home/lvotapka/Downloads/APBS-1.5-linux64/bin/apbs'
me.browndye.fhpd_numtraj = 1000

ion1 = seekr.APBS_ion('Cl-', 0.10, -1.0, 1.67)
ion2 = seekr.APBS_ion('Ca2+', 0.02, 2.0, 1.14) # define ions in the system
ion3 = seekr.APBS_ion('tris', 0.06, 1.0, 4.0)


me.browndye.apbs.ions = [ion1, ion2, ion3]
me.browndye.apbs.linear_pbe = False
me.browndye.apbs.inputgen.executable = '/home/lvotapka/Downloads/APBS-1.5-linux64/share/apbs/tools/manip/inputgen.py'

# Generate Milestones
origin = np.array([-1.536, 13.859, 16.539])
#radius_list = np.arange(2.0, 14.1, 2.0)
radius_list = [10.0, 12.0, 14.0]
vectors = [np.array([9.019, 72.142, 15.943]),]
milestones = seekr.generate_spherical_milestones(me, rec_site_atom_indices, origin, radius_list, 0, vectors, absolute=False)
print("The following milestones were created:")
seekr.print_spherical_milestone_info(milestones)
me.milestones = milestones

# Generate Filetree and Building files
seekr.generate_filetree(me)
holo_config_wet, insert_index, last_ligand_index = seekr.generate_configs(me)
ligand_heavy_indices = seekr.find_heavy_atoms(holo_config_wet, list(range(insert_index, insert_index+last_ligand_index+1)))
print("Ligand heavy indices:", ligand_heavy_indices)  # ATTN: changed

amber_settings = amber.AmberSettings()
amber_settings.leap_program = 'tleap'

amber_settings.leap_template = '''
source leaprc.protein.ff14SB
source leaprc.gaff
source leaprc.water.tip4pew
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

me.min_equil.constrained += list(range(insert_index))
me.min_equil.constrained += ligand_heavy_indices

for milestone in me.milestones:
    milestone.atom_selection_1 = me.selections.site_com_indices
    milestone.atom_selection_2 = ligand_heavy_indices
    if milestone.md:
        #amber.amber_building(me, milestone, amber_settings)

        if not milestone.openmm.prmtop_filename: continue
        # modify the file to have the correct solvent octahedron box
        parm = pmd.load_file(milestone.openmm.prmtop_filename, xyz=milestone.openmm.inpcrd_filename)

        #TODO: straighten out these units
        parm.box = np.array([64.9127105, 64.9127105, 64.9127105,109.471219, 109.471219,109.471219])
        box_vector = Quantity([[64.912710500000003, 0.0, 0.0], [-21.637568420791037, 61.200290990259163, 0.0], [-21.637568420791037, -30.600141791568205, 53.00100885481632]], unit=angstrom)
        parm.box_vectors = box_vector
        print("parm.box:", parm.box)
        print("parm.box_vectors:", parm.box_vectors)
        print("saving prmtop for milestone:", milestone.index)
        parm.save(milestone.openmm.prmtop_filename, overwrite = True)
        print("saving inpcrd for milestone:", milestone.index)
        parm.save(milestone.openmm.inpcrd_filename, overwrite = True)

        amber.create_simulation(me, milestone)
        milestone.openmm.simulation.context.setPeriodicBoxVectors([6.4912710500000003, 0.0, 0.0], [-2.1637568420791037, 6.1200290990259163, 0.0], [-2.1637568420791037, -3.0600141791568205, 5.300100885481632])

        if not me.openmm.simulation: # create a sample of a simulation file for future use
            me.openmm.simulation = milestone.openmm.simulation

#seekr.run_min_equil(me)

# save equilibration output
for milestone in me.milestones:
    if milestone.md:
        filename = amber.save_restart(me, milestone)
        milestone.openmm.umbrella_pdb_filename = filename

print("Saving all system settings for Umbrella stage and later stages.")
me.save()

print("Preparing Brownian dynamics stages...")

bd.build_bd(me)

print("Ligand heavy indices:", ligand_heavy_indices)
