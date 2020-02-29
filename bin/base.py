'''
The base module contains all the classes used in SEEKR runs

Created on May 8, 2018

@author: lvotapka
'''

import pickle as pickle
import os

class _Project():
    '''An object for generic project-level information of a SEEKR run'''
    def __init__(self):
        self.name = ''
        self.rootdir = '' # a directory to write all the MD and BD files
        self.test_mode = False # reduces the calculation time by restricting input size. Use only for debugging SEEKR
        self.md = True # whether MD is run in this calculation
        self.bd = True # whether BD is run in this calculation
        #self.k_off = False # whether files to prepare k-off calculations are generated
        self.empty_rootdir = False # if set to True will empty the contents of the rootdir when the new file tree is made

        ''' Save these for the umbrella object
        self.umbrella 100000000 # number of timesteps
        self.number_of_ens_equil_frames 10000 # number of frames to write after the ens_equil_simulations
        self.number_of_ens_equil_frames_skipped 3000
        self.extract_stride 10
        '''

class _OpenMM():
    '''An object to represent all information about the molecules. Particularly structural information.'''
    def __init__(self):
        ''' I'm making a separate program to generate the prmtop/inpcrd combo, so I'm going to assume that the user has already generated the molecule
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

class _Inputgen():
    '''Inputgen is a program for generating APBS input files.'''
    def __init__(self):
        self.executable = 'inputgen' # the location of the Inputgen executable
        self.fadd = 150 # how many angstroms around the molecule to add as a buffer before reaching the boundary in the fine grid
        self.cfac = 6 # the factor by which to create a boundary in the course grid
        self.gmemceil = 64000
        self.resolution = 0.5
        self.ionic_str = 0.15

class APBS_ion():
    '''An object to represent ions in the APBS calculation.'''
    def __init__(self, name, concentration, charge, radius):
        self.name = name # the name of the ion, can be 'Cl-', 'Ca2+', 'tris', etc.
        self.concentration = concentration # in M
        self.charge = charge # in proton charge units
        self.radius = radius # in Angstroms

class _APBS():
    '''APBS is an electrostatic Poisson-Boltzmann equation solver used in
    Browndye.'''
    def __init__(self):
        self.executable = 'apbs' # the location of the APBS executable
        self.inputgen = _Inputgen()
        self.ions = [] # a list of APBS_ion objects representing all ions in the APBS calculation
        self.linear_pbe = True # If set to True, the 'lpbe' option will be used in APBS. If set to False, the 'npbe' option will be used.

class _Browndye():
    '''An object to represent all the BrownDye molecule objects and settings.'''
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

class _Selections():
    '''An object to represent selections of ligand and receptor objects for simulation purposes.'''
    def __init__(self):
        self.lig_com_indices = [] # indices of ligand atoms to count as the center of mass
        self.site_com_indices = [] # indices to count as the center of mass of the 'site' (could be center of mass of the membrane in a permeability calculation

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

class _Min_Equil():
    '''An object for minimization and equilibration.'''
    def __init__(self):
        # Minimization variables
        self.min = True # whether to minimize the structure
        #self.min_constrained = [] # a list of atom indices
        self.min_num_steps = 0
        self.min_reporter_freq = 0 # the OpenMM 'reporter' that outputs minimization structures
        # Temperature step equilibration variables
        self.temp_equil = True # whether to run a temperature equilibration
        self.constrained = [] # list of atom indices
        self.temp_equil_temperatures = [] # a list of temperatures
        self.temp_equil_steps = 0
        self.temp_equil_reporters = []
        self.temp_equil_integrator = None # the OpenMM integrator object

class _Umbrella():
    '''An object for containing all settings related to umbrella sampling.'''
    def __init__(self):
        self.umbrella = True
        self.steps = 0
        self.energy_freq = 1
        self.traj_freq = 1
        self.force = None # the force object for umbrella sampling
        self.force_constant = 0.0
        self.integrator = None # the OpenMM integrator object
        self.reporters = [] # OpenMM reporter frequency
        self.barostat = True
        self.barostat_coeff = 25
        self.barostat_pressure = 1.0
        self.traj = []

class _Fwd_rev():
    '''An object to contain all forward-reverse stage parameters.'''
    def __init__(self):
        self.integrator = None # OpenMM integrator object
        self.reporters = [] # OpenMM reporter frequency
        self.launches_per_config = 1 # For each umbrella conformation, this represents the number of times to reinitialize the velocities and relaunch
        self.reversal_coords_pickle = ''
        self.reversal_vels_pickle = ''
        self.success_coords_pickle = ''
        self.success_vels_pickle = ''


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

    def save(self, picklename=''):
        '''Save a copy of this SEEKR calculation and all its milestone information.'''

        # need to step through the object and set all unpicklable objects to None
        self.openmm.system = None
        self.openmm.simulation = None
        self.openmm.context = None
        self.openmm.platform = None
        self.min_equil.temp_equil_integrator = None
        for milestone in self.milestones:
            milestone.openmm.system = None
            milestone.openmm.simulation = None

        if not picklename:
            picklename = os.path.join(self.project.rootdir, 'seekr_calc.pickle')
        #dill.detect.trace(True)
        #dill.detect.errors(self)
        our_file=open(picklename, 'wb')
        pickle.dump(self, our_file, protocol=-1) # protocol=-1 means to use the fastest available protocol method
        our_file.close()

    '''  # not necessary because it will be in the seekr pickle ???
    def save_milestones(self, basename='milestone.pickle'):
      for milestone in self.milestones:
        picklename = os.path.join(milestone.directory, basename)
        our_file = open(picklename, 'wb')
        pickle.dump(milestone)
    '''

def openSeekrCalc(picklename):
    our_file=open(picklename, 'rb')
    seekr_obj=pickle.load(our_file)
    our_file.close()
    return seekr_obj
