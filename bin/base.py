'''
The base module contains all the classes used in SEEKR runs

Created on May 8, 2018

@author: lvotapka
'''

import pickle
import os
import unittest, warnings
import xml.etree.ElementTree as ET
from xml.dom import minidom


from simtk.unit import kilocalorie, angstrom, mole, bar

def strBool(bool_str):
    '''
    Takes the string "true" or "false" of any case and returns a boolean object.
    '''
    if bool_str.lower() == 'true':
        return True
    elif bool_str.lower() == 'false':
        return False
    else:
        raise Exception(
            "argument for strBool must be string either 'True' or 'False'.")
    return

class _Project():
    '''An object for generic project-level information of a SEEKR run'''
    def __init__(self):
        self.name = ''
        self.rootdir = '' # a directory to write all the MD and BD files
        self.test_mode = False # reduces the calculation time by restricting 
        # input size. Use only for debugging SEEKR
        self.md = True # whether MD is run in this calculation
        self.bd = True # whether BD is run in this calculation
        #self.k_off = False # whether files to prepare k-off calculations are 
        # generated
        self.empty_rootdir = False # if set to True will empty the contents of 
        # the rootdir when the new file tree is made
        return
    
    def serialize(self, xmlProject):
        xmlProjectName = ET.SubElement(xmlProject, 'name')
        xmlProjectName.text = self.name
        xmlProjectRootdir = ET.SubElement(xmlProject, 'rootdir')
        xmlProjectRootdir.text = self.rootdir
        xmlProjectTestMode = ET.SubElement(xmlProject, 'test_mode')
        xmlProjectTestMode.text = str(self.test_mode)
        xmlProjectMd = ET.SubElement(xmlProject, 'md')
        xmlProjectMd.text = str(self.md)
        xmlProjectBd = ET.SubElement(xmlProject, 'bd')
        xmlProjectBd.text = str(self.bd)
        xmlProjectEmptyRootdir = ET.SubElement(xmlProject, 'empty_rootdir')
        xmlProjectEmptyRootdir.text = str(self.empty_rootdir)
        return
    
    def deserialize(self, xmlTree):
        self.name = xmlTree.find('name').text
        self.rootdir = xmlTree.find('rootdir').text
        self.test_mode = strBool(xmlTree.find('test_mode').text)
        self.md = strBool(xmlTree.find('md').text)
        self.bd = strBool(xmlTree.find('bd').text)
        self.empty_rootdir = strBool(xmlTree.find('empty_rootdir').text)
        return

class _OpenMM():
    '''
    An object to represent all information about the molecules. Particularly 
    structural information.
    '''
    def __init__(self):
        ''' I'm making a separate program to generate the prmtop/inpcrd combo, 
        so I'm going to assume that the user has already generated the molecule
        self.lig_pdb_filename = '' # a dry ligand PDB
        self.lig_pqr_filename = '' # the ligand PQR file
        self.rec_wet_pdb_filename = '' # A wet receptor PDB molecule
        self.rec_wet_psf_filename = '' # (only necessary for a CHARMM run)
        self.rec_dry_pdb_filename = ''
        self.rec_dry_pqr_filename = ''
        '''

        self.system = None # the OpenMM molecular system ready for MD
        self.context = None # starting context for OpenMM calculations
        self.simulation = None
        self.platform = None # OpenMM platform object
        self.properties = {} # OpenMM platform properties
        return

    def serialize(self, xmlOpenMM):
        xmlProperties = ET.SubElement(xmlOpenMM, 'properties')
        for property_key in self.properties:
            xmlProperties_i = ET.SubElement(xmlProperties, property_key)
            xmlProperties_i.text = self.properties[property_key]
        return   
       
    def deserialize(self, xmlTree):
        self.properties = {}
        for property_key in xmlTree.find('properties'):
            value = property_key.text
            self.properties[property_key.tag] = value
        return  

class _Inputgen():
    '''Inputgen is a program for generating APBS input files.'''
    def __init__(self):
        self.executable = 'inputgen' # the location of the Inputgen executable
        self.fadd = 150 # how many angstroms around the molecule to add as a 
        #buffer before reaching the boundary in the fine grid
        self.cfac = 6 # the factor to create a boundary in the course grid
        self.gmemceil = 1000
        self.resolution = 0.5
        self.ionic_str = 0.15
        return

    def serialize(self, xmlInputGen):
        xmlExecutable = ET.SubElement(xmlInputGen, 'executable')
        xmlExecutable.text = self.executable
        xmlFadd = ET.SubElement(xmlInputGen, 'fadd')
        xmlFadd.text = str(self.fadd)
        xmlCfac = ET.SubElement(xmlInputGen, 'cfac')
        xmlCfac.text = str(self.cfac)
        xmlGmemceil = ET.SubElement(xmlInputGen, 'gmemceil')
        xmlGmemceil.text = str(self.gmemceil)
        xmlResolution = ET.SubElement(xmlInputGen, 'resolution')
        xmlResolution.text = str(self.resolution)
        xmlIonicStr = ET.SubElement(xmlInputGen, 'ionic_str')
        xmlIonicStr.text = str(self.ionic_str)
        return
    
    def deserialize(self, xmlTree):
        self.executable = xmlTree.find('executable').text
        self.fadd = int(xmlTree.find('fadd').text)
        self.cfac = int(xmlTree.find('cfac').text)
        self.gmemceil = int(xmlTree.find('gmemceil').text)
        self.resolution = float(xmlTree.find('resolution').text)
        self.ionic_str = float(xmlTree.find('ionic_str').text)
        return

class APBS_ion():
    '''An object to represent ions in the APBS calculation.'''
    def __init__(self, name='', concentration=0.0, charge=0.0, radius=0.0):
        self.name = name # the name of the ion, 'Cl-', 'Ca2+', 'tris', etc.
        self.concentration = concentration # in M
        self.charge = charge # in proton charge units
        self.radius = radius # in Angstroms
        return

    def serialize(self, xmlApbs_ion):
        xmlName = ET.SubElement(xmlApbs_ion, 'name')
        xmlName.text = self.name
        xmlConcentration = ET.SubElement(xmlApbs_ion, 'concentration')
        xmlConcentration.text = str(self.concentration)
        xmlCharge = ET.SubElement(xmlApbs_ion, 'charge')
        xmlCharge.text = str(self.charge)
        xmlRadius = ET.SubElement(xmlApbs_ion, 'radius')
        xmlRadius.text = str(self.radius)
        return
    
    def deserialize(self, xmlTree):
        self.name = xmlTree.find('name').text
        self.concentration = float(xmlTree.find('concentration').text)
        self.charge = float(xmlTree.find('charge').text)
        self.radius = float(xmlTree.find('radius').text)
        return

class _APBS():
    '''APBS is an electrostatic Poisson-Boltzmann equation solver used in
    Browndye.'''
    def __init__(self):
        self.executable = 'apbs' # the location of the APBS executable
        self.inputgen = _Inputgen()
        self.ions = [] # a list of APBS_ion objects representing all ions in 
        # the APBS calculation
        self.linear_pbe = True # If set to True, the 'lpbe' option will be 
        # used in APBS. If set to False, the 'npbe' option will be used.
        return
    
    def serialize(self, xmlApbs):
        xmlExecutable = ET.SubElement(xmlApbs, 'executable')
        xmlExecutable.text = self.executable
        xmlInputGen = ET.SubElement(xmlApbs, 'inputgen')
        xmlInputGen.text = self.inputgen.serialize(xmlInputGen)
        xmlApbsIons = ET.SubElement(xmlApbs, 'apbs_ions')
        for apbs_ion in self.ions:
            xmlApbs_ion_i = ET.SubElement(xmlApbsIons, 'apbs_ion')
            apbs_ion.serialize(xmlApbs_ion_i)
        
        xmlLinear_pbe = ET.SubElement(xmlApbs, 'linear_pbe')
        xmlLinear_pbe.text = str(self.linear_pbe)
        return
    
    def deserialize(self, xmlTree):
        self.executable = xmlTree.find('executable').text
        self.inputgen.deserialize(xmlTree.find('inputgen'))
        for apbs_ion_child in xmlTree.find('apbs_ions'):
            apbs_ion = APBS_ion()
            apbs_ion.deserialize(apbs_ion_child)
            self.ions.append(apbs_ion)
        self.linear_pbe = strBool(xmlTree.find('linear_pbe').text)
        '''
        print('executable', self.executable, 
              'linear_pbe', self.linear_pbe, 
              ) 
        ''' # TODO: marked for removal
        return

class _Browndye():
    '''
    An object to represent all the BrownDye molecule objects and settings.
    '''
    def __init__(self):
        self.rec_dry_pqr_filename = ''
        self.lig_dry_pqr_filename = ''
        # BrownDye parameters
        self.b_surface_path = ''
        self.starting_lig_config = None
        self.browndye_bin = ''
        self.num_threads = 1
        self.prods_per_anchor = 1 # how many trajectories to relaunch per anchor
        # APBS parameters
        self.apbs = _APBS()
        self.ligand_is_protein = False
        self.fhpd_numtraj = 1000
        return
        
    def serialize(self, xmlBrowndye):
        xmlRec_dry_pqr_filename = ET.SubElement(xmlBrowndye, 
                                                'rec_dry_pqr_filename')
        xmlRec_dry_pqr_filename.text = self.rec_dry_pqr_filename
        xmlLig_dry_pqr_filename = ET.SubElement(xmlBrowndye, 
                                                'lig_dry_pqr_filename')
        xmlLig_dry_pqr_filename.text = self.lig_dry_pqr_filename
        xmlB_surface_path = ET.SubElement(xmlBrowndye, 'b_surface_path')
        xmlB_surface_path.text = self.b_surface_path
        xmlBrowndye_bin = ET.SubElement(xmlBrowndye, 'browndye_bin')
        xmlBrowndye_bin.text = self.browndye_bin
        xmlNum_threads = ET.SubElement(xmlBrowndye, 'num_threads')
        xmlNum_threads.text = str(self.num_threads)
        xmlProds_per_anchor = ET.SubElement(xmlBrowndye, 'prods_per_anchor')
        xmlProds_per_anchor.text = str(self.prods_per_anchor)
        xmlApbs = ET.SubElement(xmlBrowndye, 'apbs')
        xmlApbs.text = self.apbs.serialize(xmlApbs)
        xmlLigand_is_protein = ET.SubElement(xmlBrowndye, 'ligand_is_protein')
        xmlLigand_is_protein.text = str(self.ligand_is_protein)
        xmlFhpd_numtraj = ET.SubElement(xmlBrowndye, 'fhpd_numtraj')
        xmlFhpd_numtraj.text = str(self.fhpd_numtraj)
        return
    

    def deserialize(self, xmlTree):
        self.rec_dry_pqr_filename = xmlTree.find('rec_dry_pqr_filename').text
        self.lig_dry_pqr_filename = xmlTree.find('lig_dry_pqr_filename').text
        self.b_surface_path = xmlTree.find('b_surface_path').text
        self.browndye_bin = xmlTree.find('browndye_bin').text
        self.num_threads = int(xmlTree.find('num_threads').text)
        self.prods_per_anchor = int(xmlTree.find('prods_per_anchor').text)
        self.apbs.deserialize(xmlTree.find('apbs'))
        self.ligand_is_protein = strBool(xmlTree.find('ligand_is_protein').text)
        self.fhpd_numtraj = int(xmlTree.find('fhpd_numtraj').text)
        
        '''
        print('rec_dry_pqr_filename', self.rec_dry_pqr_filename, 
              'lig_dry_pqr_filename', self.lig_dry_pqr_filename, 
              'b_surface_path', self.b_surface_path, 
              'browndye_bin', self.browndye_bin, 
              'num_threads', self.num_threads, 
              'prods_per_anchor', self.prods_per_anchor, 
              'ligand_is_protein:', self.ligand_is_protein,
              'fhpd_numtraj:', self.fhpd_numtraj
              )''' # TODO: marked for removal
        return

class _Selections():
    '''
    An object to represent selections of ligand and receptor objects for 
    simulation purposes.
    '''
    def __init__(self):
        self.lig_com_indices = [] # indices of ligand atoms to count as the 
        # center of mass
        self.site_com_indices = [] # indices to count as the center of mass of 
        # the 'site' (could be center of mass of the membrane in a 
        # permeability calculation
        return
        
    def serialize(self, xmlSelections):
        xmlLig_com_indices = ET.SubElement(xmlSelections, 'lig_com_indices')
        xmlLig_com_indices.text = ', '.join(list(map(str, 
                                                     self.lig_com_indices)))
        '''
        for lig_com_index in self.lig_com_indices:
            xmlLig_com_index = ET.SubElement(xmlLig_com_indices, 
                                            'lig_com_index')
            xmlLig_com_index.text = str(lig_com_index)
        '''
        xmlSite_com_indices = ET.SubElement(xmlSelections, 'site_com_indices')
        xmlSite_com_indices.text = ', '.join(list(map(str, 
                                                      self.site_com_indices)))
        '''
        for site_com_index in self.site_com_indices:
            xmlSite_com_index = ET.SubElement(xmlSite_com_indices, 
                                                'site_com_index')
            xmlSite_com_index.text = str(site_com_index)
        '''
        return
    
    def deserialize(self, xmlTree):
        xmlLigCom_indices_str = xmlTree.find('lig_com_indices').text
        if xmlLigCom_indices_str:
            xmlLigCom_indices_str = xmlLigCom_indices_str.replace(' ', '')
            self.lig_com_indices = list(
                map(int, xmlLigCom_indices_str.split(',')))
        xmlSiteCom_indices_str = xmlTree.find('site_com_indices').text
        if xmlSiteCom_indices_str:
            xmlSiteCom_indices_str = xmlSiteCom_indices_str.replace(' ', '')
            self.site_com_indices = list(
                map(int, xmlSiteCom_indices_str.split(',')))
        '''
        print('lig_com_indices', self.lig_com_indices,
              'site_com_indices', self.site_com_indices)
        ''' # TODO: marked for removal
        return
    
class _Building():
    '''An object to represent the building process.'''
    def __init__(self):
        self.ff = ''
        self.lig_dry_pqr_filename = ''
        self.rec_wet_pdb_filename = ''
        self.rec_dry_pqr_filename = ''
        self.md_file_paths = []
        self.bd_file_paths = []
        self.config_dirlist = []
        self.reject_clashes = True
        self.ligand = None
        self.receptor_wet = None
        self.receptor_dry = None
        self.lig_resname = ''
        #self.prmtop = None
        #self.inpcrd = None
        # self.watermodel = '' ?
        # the rest needs to be filled out by the AmberPrepare script
        return

    def serialize(self, xmlBuilding):
        xmlFf = ET.SubElement(xmlBuilding, 'ff')
        xmlFf.text = self.ff
        ''' # removing these temporarily since they sometimes contain 
            # weird pdb2.py objects and they're not needed.
        xmlLig_dry_pqr_filename = ET.SubElement(xmlBuilding, 
                                                'lig_dry_pqr_filename')
        xmlLig_dry_pqr_filename.text = self.lig_dry_pqr_filename
        xmlRec_wet_pdb_filename = ET.SubElement(xmlBuilding, 
                                                'rec_wet_pdb_filename')
        xmlRec_wet_pdb_filename.text = self.rec_wet_pdb_filename
        xmlRec_dry_pqr_filename = ET.SubElement(xmlBuilding, 
                                                'rec_dry_pqr_filename')
        xmlRec_dry_pqr_filename.text = self.rec_dry_pqr_filename
        '''
        xmlMd_file_paths = ET.SubElement(xmlBuilding, 'md_file_paths')
        for md_file_path in self.md_file_paths:
            xmlMd_file_path = ET.SubElement(xmlMd_file_paths, 'md_file_path')
            xmlMd_file_path.text = md_file_path
        xmlBd_file_paths = ET.SubElement(xmlBuilding, 'bd_file_paths')
        for bd_file_path in self.bd_file_paths:
            xmlBd_file_path = ET.SubElement(xmlBd_file_paths, 'bd_file_path')
            xmlBd_file_path.text = bd_file_path
        xmlConfig_dirlist = ET.SubElement(xmlBuilding, 'config_dirlist')
        for config_dir in self.config_dirlist:
            xmlConfig_dir = ET.SubElement(xmlConfig_dirlist, 'config_dir')
            xmlConfig_dir.text = config_dir
        xmlReject_clashes = ET.SubElement(xmlBuilding, 'reject_clashes')
        xmlReject_clashes.text = str(self.reject_clashes)
        xmlLig_resname = ET.SubElement(xmlBuilding, 'lig_resname')
        xmlLig_resname.text = self.lig_resname
        return
    
    def deserialize(self, xmlTree):
        self.ff = xmlTree.find('ff').text
        self.md_file_paths = []
        for md_file_path in xmlTree.find('md_file_paths'):
            self.md_file_paths.append(md_file_path.text)
        self.bd_file_paths = []
        for bd_file_path in xmlTree.find('bd_file_paths'):
            self.bd_file_paths.append(bd_file_path.text)
        self.config_dirlist = []
        for config_dir in xmlTree.find('config_dirlist'):
            self.config_dirlist.append(config_dir.text)
        self.reject_clashes = strBool(xmlTree.find('reject_clashes').text)
        self.lig_resname = xmlTree.find('lig_resname').text
        '''
        print('ff', self.ff,
              'md_file_paths', self.md_file_paths,
              'bd_file_paths', self.bd_file_paths,
              'config_dirlist', self.config_dirlist,
              'reject_clashes', self.reject_clashes,
              'lig_resname', self.lig_resname
              )
        ''' # TODO: marked for removal
        return

class _Min_Equil():
    '''An object for minimization and equilibration.'''
    def __init__(self):
        # Minimization variables
        self.min = True # whether to minimize the structure
        #self.min_constrained = [] # a list of atom indices
        self.min_num_steps = 0
        self.min_reporter_freq = 0 # the OpenMM 'reporter' that outputs 
        # minimization structures
        # Temperature step equilibration variables
        self.temp_equil = True # whether to run a temperature equilibration
        self.constrained = [] # list of atom indices
        self.temp_equil_temperatures = [] # a list of temperatures
        self.temp_equil_steps = 0
        self.temp_equil_reporters = []
        self.temp_equil_integrator = None # the OpenMM integrator object
        return
    
    def serialize(self, xmlMin_equil):
        xmlMin = ET.SubElement(xmlMin_equil, 'min')
        xmlMin.text = str(self.min)
        xmlMin_num_steps = ET.SubElement(xmlMin_equil, 'min_num_steps')
        xmlMin_num_steps.text = str(self.min_num_steps)
        xmlMin_reporter_freq = ET.SubElement(xmlMin_equil, 'min_reporter_freq')
        xmlMin_reporter_freq.text = str(self.min_reporter_freq)
        xmlTemp_equil = ET.SubElement(xmlMin_equil, 'temp_equil')
        xmlTemp_equil.text = str(self.temp_equil)
        xmlConstrained = ET.SubElement(xmlMin_equil, 'constrained')
        xmlConstrained.text = ', '.join(list(map(str, self.constrained)))
        '''
        for constrained_index in self.constrained:
            xmlConstrained_index = ET.SubElement(xmlConstrained, 
                                                'constrained_index')
            xmlConstrained_index.text = str(constrained_index)
        '''
        xmlTemp_equil_temperatures = ET.SubElement(xmlMin_equil, 
                                                   'temp_equil_temperatures')
        xmlTemp_equil_temperatures.text = ', '.join(list(map(str, 
            self.temp_equil_temperatures)))
        '''
        for temp_equil_temperature in self.temp_equil_temperatures:
            xmlTemp_equil_temperature = ET.SubElement(
                xmlTemp_equil_temperatures, 'temperature')
            xmlTemp_equil_temperature.text = str(temp_equil_temperature)
        '''
        xmlTemp_equil_steps = ET.SubElement(xmlMin_equil, 'temp_equil_steps')
        xmlTemp_equil_steps.text = str(self.temp_equil_steps)
        return
    
    def deserialize(self, xmlTree):
        self.min = strBool(xmlTree.find('min').text)
        self.min_num_steps = int(xmlTree.find('min_num_steps').text)
        self.min_reporter_freq = int(xmlTree.find('min_reporter_freq').text)
        self.temp_equil = strBool(xmlTree.find('temp_equil').text)
        xmlConstrained_str = xmlTree.find('constrained').text
        if xmlConstrained_str:
            xmlConstrained_str = xmlConstrained_str.replace(' ', '')
            self.constrained = list(
                map(int, xmlConstrained_str.split(',')))
            
        xmlTemp_equil_temperatures_str = \
            xmlTree.find('temp_equil_temperatures').text
        if xmlTemp_equil_temperatures_str:
            xmlTemp_equil_temperatures_str = \
                xmlTemp_equil_temperatures_str.replace(' ', '')
            self.temp_equil_temperatures = list(
                map(float, xmlTemp_equil_temperatures_str.split(',')))
        
        self.temp_equil_steps = int(xmlTree.find('temp_equil_steps').text)
        return
    
class _Umbrella():
    '''An object for containing all settings related to umbrella sampling.'''
    def __init__(self):
        self.umbrella = True
        self.steps = 0
        self.energy_freq = 1
        self.traj_freq = 1
        self.force = None # the force object for umbrella sampling
        self.force_constant = None
        self.integrator = None # the OpenMM integrator object
        self.reporters = [] # OpenMM reporter frequency
        self.barostat = True
        self.barostat_coeff = 0
        self.barostat_pressure = None
        self.traj = []
        return
    
    def serialize(self, xmlUmbrella):
        xmlRunUmbrella = ET.SubElement(xmlUmbrella, 'umbrella')
        xmlRunUmbrella.text = str(self.umbrella)
        xmlSteps = ET.SubElement(xmlUmbrella, 'steps')
        xmlSteps.text = str(self.steps)
        xmlEnergy_freq = ET.SubElement(xmlUmbrella, 'energy_freq')
        xmlEnergy_freq.text = str(self.energy_freq)
        xmlTraj_freq = ET.SubElement(xmlUmbrella, 'traj_freq')
        xmlTraj_freq.text = str(self.traj_freq)
        xmlForce_constant = ET.SubElement(xmlUmbrella, 'force_constant')
        if self.force_constant is not None:
            xmlForce_constant.text = str(self.force_constant.value_in_unit(
                kilocalorie/(angstrom**2*mole)))
        else:
            xmlForce_constant.text = None
        xmlBarostat = ET.SubElement(xmlUmbrella, 'barostat')
        xmlBarostat.text = str(self.barostat)
        xmlBarostat_coeff = ET.SubElement(xmlUmbrella, 'barostat_coeff')
        xmlBarostat_coeff.text = str(self.barostat_coeff)
        xmlBarostat_pressure = ET.SubElement(xmlUmbrella, 'barostat_pressure')
        if self.barostat_pressure is not None:
            xmlBarostat_pressure.text = str(
                self.barostat_pressure.value_in_unit(bar))
        else:
            xmlBarostat_pressure.text = None
        return
    
    def deserialize(self, xmlTree):
        self.umbrella = strBool(xmlTree.find('umbrella').text)
        self.steps = int(xmlTree.find('steps').text)
        self.energy_freq = int(xmlTree.find('energy_freq').text)
        self.traj_freq = int(xmlTree.find('traj_freq').text)
        if xmlTree.find('force_constant').text is not None:
            force_constant_val = float(xmlTree.find('force_constant').text)
            self.force_constant = force_constant_val * \
                kilocalorie/(angstrom**2*mole)
        else:
            self.force_constant = None
        self.barostat = strBool(xmlTree.find('barostat').text)
        self.barostat_coeff = float(xmlTree.find('barostat_coeff').text)
        if xmlTree.find('barostat_pressure').text is not None:
            barostat_pressure_val = float(xmlTree.find('barostat_pressure').text)
            self.barostat_pressure = barostat_pressure_val * bar
        else:
            self.barostat_pressure = None
        return

class _Fwd_rev():
    '''An object to contain all forward-reverse stage parameters.'''
    def __init__(self):
        self.integrator = None # OpenMM integrator object
        self.reporters = [] # OpenMM reporter frequency
        self.launches_per_config = 1 # For each umbrella conformation, this 
        # represents the number of times to reinitialize the velocities and relaunch
        self.reversal_coords_pickle = ''
        self.reversal_vels_pickle = ''
        self.success_coords_pickle = ''
        self.success_vels_pickle = ''
        return
    
    def serialize(self, xmlFwd_rev):
        xmlLaunches_per_config = ET.SubElement(xmlFwd_rev, 
                                               'launches_per_config')
        xmlLaunches_per_config.text = str(self.launches_per_config)
        xmlReversal_coords_pickle = ET.SubElement(xmlFwd_rev, 
                                                  'reversal_coords_pickle')
        xmlReversal_coords_pickle.text = str(self.reversal_coords_pickle)
        xmlReversal_vels_pickle = ET.SubElement(xmlFwd_rev, 
                                                'reversal_vels_pickle')
        xmlReversal_vels_pickle.text = str(self.reversal_vels_pickle)
        xmlSuccess_coords_pickle = ET.SubElement(xmlFwd_rev, 
                                                 'success_coords_pickle')
        xmlSuccess_coords_pickle.text = str(self.success_coords_pickle)
        xmlSuccess_vels_pickle = ET.SubElement(xmlFwd_rev, 
                                               'success_vels_pickle')
        xmlSuccess_vels_pickle.text = str(self.success_vels_pickle)
        return
    
    def deserialize(self, xmlTree):
        self.launches_per_config = int(xmlTree.find('launches_per_config').text)
        self.reversal_coords_pickle = xmlTree.find(
            'reversal_coords_pickle').text
        self.reversal_vels_pickle = xmlTree.find(
            'reversal_vels_pickle').text
        self.success_coords_pickle = xmlTree.find(
            'success_coords_pickle').text
        self.success_vels_pickle = xmlTree.find(
            'success_vels_pickle').text
        return

class SeekrCalculation():
    '''An encapsulating class for all settings within a SEEKR calculation.'''
    def __init__(self):
        '''Initialize the SEEKR object'''
        self.master_temperature = 0
        self.milestones = [] # a list of milestone objects
        self.project = _Project() # project object
        self.openmm = _OpenMM() # OpenMM settings object
        self.browndye = _Browndye() # Browndye settings object
        self.selections = _Selections() # selections settings
        self.building = _Building()
        self.min_equil = _Min_Equil()
        self.umbrella_stage = _Umbrella()
        self.fwd_rev_stage = _Fwd_rev()
        return
        
    def deserialize(self, root):
        import milestones
        self.master_temperature = float(root.find('master_temperature').text)
        self.milestones = milestones.deserialize_milestones(
            root.find('milestones'))
        self.project.deserialize(root.find('project'))
        self.openmm.deserialize(root.find('openmm'))
        self.browndye.deserialize(root.find('browndye'))
        self.selections.deserialize(root.find('selections'))
        self.building.deserialize(root.find('building'))
        self.min_equil.deserialize(root.find('min_equil'))
        self.umbrella_stage.deserialize(root.find('umbrella_stage'))
        self.fwd_rev_stage.deserialize(root.find('fwd_rev_stage'))
        return

    def save(self, picklename=''):
        '''
        Save a copy of this SEEKR calculation and all its milestone 
        information.
        '''

        root = ET.Element('seekr_calc')
        xmlMasterTemp = ET.SubElement(root, 'master_temperature')
        xmlMasterTemp.text = str(self.master_temperature)
        xmlProject = ET.SubElement(root, 'project')
        self.project.serialize(xmlProject)
        
        xmlMilestones = ET.SubElement(root, 'milestones')
        for milestone in self.milestones:
            xmlMilestone = ET.SubElement(xmlMilestones, 'milestone')
            milestone.serialize(xmlMilestone)
        
        xmlOpenMM = ET.SubElement(root, 'openmm')
        self.openmm.serialize(xmlOpenMM)
        xmlBrowndye = ET.SubElement(root, 'browndye')
        self.browndye.serialize(xmlBrowndye)
        xmlSelections = ET.SubElement(root, 'selections')
        self.selections.serialize(xmlSelections)
        xmlBuilding = ET.SubElement(root, 'building')
        self.building.serialize(xmlBuilding)
        xmlMinEquil = ET.SubElement(root, 'min_equil')
        self.min_equil.serialize(xmlMinEquil)
        xmlUmbrella = ET.SubElement(root, 'umbrella_stage')
        self.umbrella_stage.serialize(xmlUmbrella)
        xmlFwd_rev = ET.SubElement(root, 'fwd_rev_stage')
        self.fwd_rev_stage.serialize(xmlFwd_rev)
        xmlstr = minidom.parseString(ET.tostring(root)).toprettyxml(
            indent="   ")
        

        if not picklename:
            picklename = os.path.join(self.project.rootdir, 'seekr_calc.xml')
        #dill.detect.trace(True)
        #dill.detect.errors(self)
        our_file=open(picklename, 'w')
        our_file.write(xmlstr)
        our_file.close()
        return

    '''  # not necessary because it will be in the seekr pickle ???
    def save_milestones(self, basename='milestone.pickle'):
      for milestone in self.milestones:
        picklename = os.path.join(milestone.directory, basename)
        our_file = open(picklename, 'wb')
        pickle.dump(milestone)
    '''

def openSeekrCalcPickle(picklename):
    warnings.warn('Pickles are being deprecated for OpenSeekr. '\
                             'Please convert SeekrCalculation pickle files to '\
                             'XML using convert_pickle.py', FutureWarning)
    our_file=open(picklename, 'rb')
    seekr_obj=pickle.load(our_file)
    our_file.close()
    return seekr_obj

def openSeekrCalc(xmlFileName):
    tree = ET.parse(xmlFileName)
    root = tree.getroot()
    seekrCalc = SeekrCalculation()
    seekrCalc.deserialize(root)
    return seekrCalc

class Test_base(unittest.TestCase):
    # several test cases to ensure the functions in this module are working
    def test_pickle(self): # test this function
        filename = '/tmp/testdump.xml'
        myobj = SeekrCalculation()
        myobj.save(filename)
        tmp_testfile = '/home/lvotapka/tryp_testXML/seekr_calc.xml' #TODO: marked for removal
        newobj = openSeekrCalc(tmp_testfile)
        #newobj = openSeekrCalc(filename)
        return
    
    def test_strBool(self):
        self.assertEqual(strBool('True'), True)
        self.assertEqual(strBool('true'), True)
        self.assertEqual(strBool('TRUE'), True)
        self.assertEqual(strBool('False'), False)
        self.assertEqual(strBool('false'), False)
        self.assertEqual(strBool('FALSE'), False)
        with self.assertRaises(Exception):
            strBool('balderdash')
    
if __name__ == "__main__":
    unittest.main()
