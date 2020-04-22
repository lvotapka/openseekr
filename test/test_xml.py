'''
Created on Apr 7, 2020

Tests the serialize and deserialize functions in the SeekrCalculation objects
and the Milestone objects.

@author: lvotapka
'''

import base, milestones
from simtk.unit import Quantity, kilocalories_per_mole, angstroms, nanometer
from simtk.unit import bar
import numpy as np

def assertEq(value1, value2):
    if value1 != value2:
        print("ERROR - values unequal:", value1, '!=', value2)
    return

def assertEqQuantities(quantity1, quantity2, TOL=1e-8):
    quantityUnit = quantity1.unit
    newQuantity = quantity1 - quantity2
    tolQuantity = Quantity(TOL, quantityUnit)
    if newQuantity**2 > tolQuantity**2:
        print("ERROR - values unequal:", quantity1, '!=', quantity2)
    return

def make_seekr_calculation():
    seekr = base.SeekrCalculation()
    seekr.master_temperature = 312.7
    # Project
    seekr.project.name = 'testProjectName'
    seekr.project.rootdir = '/path/to/test/rootdir'
    seekr.project.test_mode = True 
    seekr.project.md = True
    seekr.project.bd = True
    seekr.project.empty_rootdir = False
    # OpenMM
    seekr.openmm.properties = {'CudaDeviceIndex':'1', 'CudaPrecision':'mixed'}
    # Browndye
    seekr.browndye.rec_dry_pqr_filename = '/path/to/test_rec_dry_pqr_filename'
    seekr.browndye.lig_dry_pqr_filename = '/path/to/test_lig_dry_pqr_filename'
    seekr.browndye.b_surface_path = '/path/to/test_b_surface_path'
    seekr.browndye.browndye_bin = '/path/to/test/browndye_bin'
    seekr.browndye.num_threads = 6
    seekr.browndye.prods_per_anchor = 100
    seekr.browndye.ligand_is_protein = False
    seekr.browndye.fhpd_numtraj = 16325
    # Browndye.apbs
    seekr.browndye.apbs.executable = '/path/to/test_apbs_executable'
    seekr.browndye.apbs.ions = [ 
        base.APBS_ion('Cl-', 0.15, -1.0, 1.98),
        base.APBS_ion('Na+', 0.15, 1.0, 1.16)
        ]
    seekr.browndye.apbs.linear_pbe = True 
    # Browndye.apbs.inputgen
    seekr.browndye.apbs.inputgen.executable = '/path/to/test_inputgen' 
    seekr.browndye.apbs.inputgen.fadd = 257 
    seekr.browndye.apbs.inputgen.cfac = 652 
    seekr.browndye.apbs.inputgen.gmemceil = 32000
    seekr.browndye.apbs.inputgen.resolution = 0.75
    seekr.browndye.apbs.inputgen.ionic_str = 0.25
    # Selections
    seekr.selections.lig_com_indices = [3221, 3222, 3223, 3224, 3225, 3226, 
                                        3227, 3228, 3229]
    seekr.selections.site_com_indices = [2478, 2489, 2499, 2535, 2718, 2745, 
                                         2769, 2787, 2794, 2867, 2926]
    # Building
    seekr.building.ff = 'testFF'
    seekr.building.lig_dry_pqr_filename = \
        '/path/to/test_building_lig_dry_pqr_filename'
    seekr.building.rec_wet_pdb_filename = \
        '/path/to/test_building_rec_wet_pdb_filename'
    seekr.building.rec_dry_pqr_filename = \
        '/path/to/test_building_rec_dry_pqr_filename'
    seekr.building.md_file_paths = ['test_md_path1', 'test_md_path2']
    seekr.building.bd_file_paths = ['test_bd_path1', 'test_bd_path2']
    seekr.building.config_dirlist = ['test_config_dir1', 'test_config_dir2',
                                     'test_config_dir3']
    seekr.building.reject_clashes = True
    seekr.building.lig_resname = 'test_lig_resname'
    # Min_equil
    seekr.min_equil.min = True
    seekr.min_equil.min_num_steps = 1234
    seekr.min_equil.min_reporter_freq = 34
    seekr.min_equil.temp_equil = True
    seekr.min_equil.constrained = [123, 456, 789, 12, 345]
    seekr.min_equil.temp_equil_temperatures = [300, 315, 330, 315, 300]
    seekr.min_equil.temp_equil_steps = 234
    # Umbrella
    seekr.umbrella_stage.umbrella = True
    seekr.umbrella_stage.steps = 100234
    seekr.umbrella_stage.energy_freq = 56
    seekr.umbrella_stage.traj_freq = 34
    seekr.umbrella_stage.force_constant = \
        39.0*kilocalories_per_mole/angstroms**2
    seekr.umbrella_stage.barostat = True
    seekr.umbrella_stage.barostat_coeff = 98
    seekr.umbrella_stage.barostat_pressure = 1.01*bar
    # Fwd_rev
    seekr.fwd_rev_stage.launches_per_config = 11
    seekr.fwd_rev_stage.reversal_coords_pickle = '/path/to/test_reversal_coords'
    seekr.fwd_rev_stage.reversal_vels_pickle = '/path/to/test_reversal_vels'
    seekr.fwd_rev_stage.success_coords_pickle = '/path/to/test_success_coords'
    seekr.fwd_rev_stage.success_vels_pickle = '/path/to/test_success_vels'
    
    for i in range(3):
        milestone = milestones.Concentric_Spherical_Milestone(i, 0)
        milestone.fullname = 'test_milestone_name_%d' % i
        milestone.directory = 'test_milestone_directory_%d' % i
        milestone.anchor = np.array([12.3, 4.56, 78.9])
        milestone.center_atom_indices = [123, 234, 345, 456] 
        milestone.center_vec = [12.3, 4.56, 78.9]
        milestone.radius = 2.3 * i
        milestone.wet_holo_filename = '/path/to/test_wet_holo_filename'
        milestone.dry_holo_filename = '/path/to/test_dry_holo_filename'
        milestone.atom_selection_1 = [56, 78, 89]
        milestone.atom_selection_2 = [34, 76, 87]
        milestone.building_box_vectors = Quantity([[64.913, 0.0, 0.0], 
                                          [-21.63756, 61.20029, 0.0], 
                                          [-21.63756, -30.6001417, 53.001008]], 
                                          unit=angstroms)
        milestone.min_equil_box_vectors = Quantity([[64.913, 0.0, 0.0], 
                                          [-20.63756, 61.20029, 0.0], 
                                          [-21.63756, -30.6001417, 53.001008]], 
                                          unit=angstroms)
        milestone.umbrella_box_vectors = Quantity([[64.413, 0.0, 0.0], 
                                          [-21.63756, 61.20029, 0.0], 
                                          [-21.63756, -30.6001417, 53.001008]], 
                                          unit=angstroms)
        milestone.fwd_rev_box_vectors = Quantity([[64.913, 0.0, 0.0], 
                                          [-21.63756, 61.20029, 0.0], 
                                          [-21.54756, -30.6001417, 53.001008]], 
                                          unit=angstroms)
        # Milestone_system
        milestone.openmm = milestones.Milestone_System()
        milestone.openmm.wet_holo_pdb_filename = \
            '/path/to/test_wet_holo_pdb_filename'
        milestone.openmm.dry_holo_pdb_filename = \
            '/path/to/test_dry_holo_pdb_filename'
        milestone.openmm.prmtop_filename = \
            '/path/to/test_prmtop_filename'
        milestone.openmm.inpcrd_filename = \
            '/path/to/test_inpcrd_filename'
        milestone.openmm.umbrella_pdb_filename = \
            '/path/to/test_umbrella_pdb_filename'
        seekr.milestones.append(milestone)
    
    for i in range(2):
        seekr.milestones[i+1].neighbors.append(seekr.milestones[i])
        
    for i in range(2):
        seekr.milestones[i].neighbors.append(seekr.milestones[i+1])
        
    return seekr
    
def verify_seekr_identical(seekr1, seekr2):
    assertEq(seekr1.master_temperature, seekr2.master_temperature)
    # Project
    assertEq(seekr1.project.name, seekr2.project.name)
    assertEq(seekr1.project.rootdir, seekr2.project.rootdir)
    assertEq(seekr1.project.test_mode, seekr2.project.test_mode)
    assertEq(seekr1.project.md, seekr2.project.md)
    assertEq(seekr1.project.bd, seekr2.project.bd)
    assertEq(seekr1.project.empty_rootdir, seekr2.project.empty_rootdir)
    # OpenMM
    assertEq(seekr1.openmm.properties, seekr2.openmm.properties)
    # Browndye
    assertEq(seekr1.browndye.rec_dry_pqr_filename, 
             seekr2.browndye.rec_dry_pqr_filename)
    assertEq(seekr1.browndye.lig_dry_pqr_filename, 
             seekr2.browndye.lig_dry_pqr_filename)
    assertEq(seekr1.browndye.b_surface_path, seekr2.browndye.b_surface_path)
    assertEq(seekr1.browndye.browndye_bin, seekr2.browndye.browndye_bin)
    assertEq(seekr1.browndye.num_threads, seekr2.browndye.num_threads)
    assertEq(seekr1.browndye.prods_per_anchor, 
             seekr2.browndye.prods_per_anchor)
    assertEq(seekr1.browndye.ligand_is_protein, 
             seekr2.browndye.ligand_is_protein)
    assertEq(seekr1.browndye.fhpd_numtraj, seekr2.browndye.fhpd_numtraj)
    # Browndye.apbs
    assertEq(seekr1.browndye.apbs.executable, seekr2.browndye.apbs.executable)
    assertEq(len(seekr1.browndye.apbs.ions), len(seekr2.browndye.apbs.ions))
    for i in range(len(seekr1.browndye.apbs.ions)):
        assertEq(seekr1.browndye.apbs.ions[i].name, 
                 seekr2.browndye.apbs.ions[i].name)
        assertEq(seekr1.browndye.apbs.ions[i].concentration, 
                 seekr2.browndye.apbs.ions[i].concentration)
        assertEq(seekr1.browndye.apbs.ions[i].charge, 
                 seekr2.browndye.apbs.ions[i].charge)
        assertEq(seekr1.browndye.apbs.ions[i].radius, 
                 seekr2.browndye.apbs.ions[i].radius)
    # Browndye.apbs.inputgen
    assertEq(seekr1.browndye.apbs.inputgen.executable, 
             seekr2.browndye.apbs.inputgen.executable)
    assertEq(seekr1.browndye.apbs.inputgen.fadd, 
             seekr2.browndye.apbs.inputgen.fadd)
    assertEq(seekr1.browndye.apbs.inputgen.cfac, 
             seekr2.browndye.apbs.inputgen.cfac)
    assertEq(seekr1.browndye.apbs.inputgen.gmemceil, 
             seekr2.browndye.apbs.inputgen.gmemceil)
    assertEq(seekr1.browndye.apbs.inputgen.resolution, 
             seekr2.browndye.apbs.inputgen.resolution)
    assertEq(seekr1.browndye.apbs.inputgen.ionic_str, 
             seekr2.browndye.apbs.inputgen.ionic_str)
    # Selections
    assertEq(seekr1.selections.lig_com_indices, 
             seekr2.selections.lig_com_indices)
    assertEq(seekr1.selections.site_com_indices, 
             seekr2.selections.site_com_indices)
    # Building
    assertEq(seekr1.building.ff, seekr2.building.ff)
    #assertEq(seekr1.building.lig_dry_pqr_filename,
    #         seekr2.building.lig_dry_pqr_filename)
    #assertEq(seekr1.building.rec_wet_pdb_filename, 
    #         seekr2.building.rec_wet_pdb_filename)
    #assertEq(seekr1.building.rec_dry_pqr_filename, 
    #         seekr2.building.rec_dry_pqr_filename)
    assertEq(len(seekr1.building.md_file_paths), 
             len(seekr2.building.md_file_paths))
    for i in range(len(seekr1.building.md_file_paths)):
        assertEq(seekr1.building.md_file_paths[i], 
                 seekr2.building.md_file_paths[i])
    assertEq(len(seekr1.building.bd_file_paths), 
             len(seekr2.building.bd_file_paths))
    for i in range(len(seekr1.building.bd_file_paths)):
        assertEq(seekr1.building.bd_file_paths[i], 
                 seekr2.building.bd_file_paths[i])
    assertEq(len(seekr1.building.config_dirlist), 
             len(seekr2.building.config_dirlist))
    for i in range(len(seekr1.building.config_dirlist)):
        assertEq(seekr1.building.config_dirlist[i], 
                 seekr2.building.config_dirlist[i])
    assertEq(seekr1.building.reject_clashes, seekr2.building.reject_clashes)
    assertEq(seekr1.building.lig_resname, seekr2.building.lig_resname)
    # Min_equil
    assertEq(seekr1.min_equil.min, seekr2.min_equil.min)
    assertEq(seekr1.min_equil.min_num_steps, seekr2.min_equil.min_num_steps)
    assertEq(seekr1.min_equil.min_reporter_freq, 
             seekr2.min_equil.min_reporter_freq)
    assertEq(seekr1.min_equil.temp_equil, seekr2.min_equil.temp_equil)
    assertEq(len(seekr1.min_equil.constrained), 
             len(seekr2.min_equil.constrained))
    for i in range(len(seekr1.min_equil.constrained)):
        assertEq(seekr1.min_equil.constrained[i], 
                 seekr2.min_equil.constrained[i])
    assertEq(len(seekr1.min_equil.temp_equil_temperatures), 
             len(seekr2.min_equil.temp_equil_temperatures))
    for i in range(len(seekr1.min_equil.temp_equil_temperatures)):
        assertEq(seekr1.min_equil.temp_equil_temperatures[i], 
                 seekr2.min_equil.temp_equil_temperatures[i])
    assertEq(seekr1.min_equil.temp_equil_steps, 
             seekr2.min_equil.temp_equil_steps)
    # Umbrella
    assertEq(seekr1.umbrella_stage.umbrella, 
             seekr2.umbrella_stage.umbrella)
    assertEq(seekr1.umbrella_stage.steps, 
             seekr2.umbrella_stage.steps)
    assertEq(seekr1.umbrella_stage.energy_freq, 
             seekr2.umbrella_stage.energy_freq)
    assertEq(seekr1.umbrella_stage.traj_freq, 
             seekr2.umbrella_stage.traj_freq)
    assertEqQuantities(seekr1.umbrella_stage.force_constant, 
                       seekr2.umbrella_stage.force_constant)
    assertEq(seekr1.umbrella_stage.barostat, 
             seekr2.umbrella_stage.barostat)
    assertEq(seekr1.umbrella_stage.barostat_coeff, 
             seekr2.umbrella_stage.barostat_coeff)
    assertEqQuantities(seekr1.umbrella_stage.barostat_pressure, 
                       seekr2.umbrella_stage.barostat_pressure)
    # Fwd_rev
    assertEq(seekr1.fwd_rev_stage.launches_per_config, 
             seekr2.fwd_rev_stage.launches_per_config)
    assertEq(seekr1.fwd_rev_stage.reversal_coords_pickle, 
             seekr2.fwd_rev_stage.reversal_coords_pickle)
    assertEq(seekr1.fwd_rev_stage.reversal_vels_pickle, 
             seekr2.fwd_rev_stage.reversal_vels_pickle)
    assertEq(seekr1.fwd_rev_stage.success_coords_pickle, 
             seekr2.fwd_rev_stage.success_coords_pickle)
    assertEq(seekr1.fwd_rev_stage.success_vels_pickle, 
             seekr2.fwd_rev_stage.success_vels_pickle)
    
    # milestones
    assertEq(len(seekr1.milestones), len(seekr2.milestones))
    for i in range(len(seekr1.milestones)):
        assertEq(seekr1.milestones[i].fullname, seekr2.milestones[i].fullname)
        assertEq(seekr1.milestones[i].directory, seekr2.milestones[i].directory)
        for j in range(3):
            assertEq(seekr1.milestones[i].anchor[j], \
                     seekr2.milestones[i].anchor[j])
        assertEq(seekr1.milestones[i].center_atom_indices, 
                 seekr2.milestones[i].center_atom_indices)
        assertEq(seekr1.milestones[i].center_vec, \
                 seekr2.milestones[i].center_vec)
        assertEq(seekr1.milestones[i].radius, seekr2.milestones[i].radius)
        assertEq(seekr1.milestones[i].wet_holo_filename, 
                 seekr2.milestones[i].wet_holo_filename)
        assertEq(seekr1.milestones[i].dry_holo_filename, 
                 seekr2.milestones[i].dry_holo_filename)
        assertEq(seekr1.milestones[i].atom_selection_1, 
                 seekr2.milestones[i].atom_selection_1)
        assertEq(seekr1.milestones[i].atom_selection_2, 
                 seekr2.milestones[i].atom_selection_2)
        
        if (seekr1.milestones[i].building_box_vectors is not None) and \
                (seekr2.milestones[i].building_box_vectors is not None):
            for j in range(3):
                for k in range(3):
                    assertEqQuantities(seekr1.milestones[i].building_box_vectors[j][k], 
                               seekr2.milestones[i].building_box_vectors[j][k])
        else:
            assert seekr1.milestones[i].building_box_vectors is None
            assert seekr2.milestones[i].building_box_vectors is None
        
        if (seekr1.milestones[i].min_equil_box_vectors is not None) and \
                (seekr2.milestones[i].min_equil_box_vectors is not None):
            for j in range(3):
                for k in range(3):
                    assertEqQuantities(seekr1.milestones[i].min_equil_box_vectors[j][k], 
                               seekr2.milestones[i].min_equil_box_vectors[j][k])
        else:
            assert seekr1.milestones[i].min_equil_box_vectors is None
            assert seekr2.milestones[i].min_equil_box_vectors is None
                
        if seekr1.milestones[i].umbrella_box_vectors is not None:
            for j in range(3):
                for k in range(3):
                    assertEqQuantities(seekr1.milestones[i].umbrella_box_vectors[j][k], 
                               seekr2.milestones[i].umbrella_box_vectors[j][k])
        else:
            assert seekr1.milestones[i].umbrella_box_vectors is None
            assert seekr2.milestones[i].umbrella_box_vectors is None
                
        if seekr1.milestones[i].fwd_rev_box_vectors is not None:
            for j in range(3):
                for k in range(3):
                    assertEqQuantities(seekr1.milestones[i].fwd_rev_box_vectors[j][k], 
                               seekr2.milestones[i].fwd_rev_box_vectors[j][k])
        else:
            assert seekr1.milestones[i].fwd_rev_box_vectors is None
            assert seekr2.milestones[i].fwd_rev_box_vectors is None
        
        # more here
        # Milestone_system
        assertEq(seekr1.milestones[i].openmm.wet_holo_pdb_filename, 
                 seekr2.milestones[i].openmm.wet_holo_pdb_filename)
        assertEq(seekr1.milestones[i].openmm.dry_holo_pdb_filename, 
                 seekr2.milestones[i].openmm.dry_holo_pdb_filename)
        assertEq(seekr1.milestones[i].openmm.prmtop_filename, 
                 seekr2.milestones[i].openmm.prmtop_filename)
        assertEq(seekr1.milestones[i].openmm.inpcrd_filename, 
                 seekr2.milestones[i].openmm.inpcrd_filename)
        assertEq(seekr1.milestones[i].openmm.umbrella_pdb_filename, 
                 seekr2.milestones[i].openmm.umbrella_pdb_filename)
        
    assertEq(seekr1.milestones[0].neighbors[0].index,1)
    assertEq(seekr1.milestones[-1].neighbors[0].index, 
             len(seekr1.milestones)-2)
    for i in range(1, len(seekr1.milestones)-1):
        assertEq(seekr1.milestones[i].neighbors[0].index, i-1)
        assertEq(seekr1.milestones[i].neighbors[1].index, i+1)
    
    return


if __name__ == "__main__":
    xmlFilename = './seekrTest.xml'
    seekr1 = make_seekr_calculation()
    seekr1.save(xmlFilename)
    seekr2 = base.openSeekrCalc(xmlFilename)
    verify_seekr_identical(seekr1, seekr2)
    