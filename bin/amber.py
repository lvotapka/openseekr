'''
Merge ligand and receptor files to write parameter/topology/coordinate files
for Amber simulations.

Created on May 10, 2018

@author: lvotapka
@contributor Ilker Deveci
'''

import os, string
from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from copy import deepcopy
import parmed

verbose = True

class AmberSettings():
    def __init__(self):
        self.leap_template = ""
        self.leap_program = 'tleap'

def amber_building(seekrcalc, milestone, amber_settings):
    '''the pre-minimization procedure for amber ff simulations
    Input:
     - seekrcalc: The SeekrCalculation object that contains all the settings for
         the SEEKR calculation.
     - milestone: the Milestone() object to prepare the AMBER system for
     - amber_settings: the AmberSettings() object that defines settings for the
         tleap calculation.
    Output:
     - None
    '''
    i = milestone.index
    building = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'building')
    working_pdb_base = 'holo' #settings['working_pdb_base']
    pdbfile = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'holo_wet.pdb')
    leap_program = amber_settings.leap_program
    if not os.path.exists(pdbfile):
        print("Skipping milestone %d: no holo structure detected at %s." % (i, pdbfile))
        return
    prmtop = os.path.join(building,working_pdb_base+'.parm7')
    inpcrd = os.path.join(building,working_pdb_base+'.rst7')
    newpdb = os.path.join(building,working_pdb_base+'_leap.pdb')

    leap_params = { 'HOLO_WET_PDB':pdbfile,
                    'PRMTOP':prmtop,
                    'INPCRD':inpcrd,
                    'LEAP_OUTPUT_PDB':newpdb,
      }
    leap_template = string.Template(amber_settings.leap_template)
    leap_string = leap_template.substitute(leap_params)

    leapfilename = os.path.join(building,'anchor.leap')

    if os.path.exists(prmtop) and os.path.exists(inpcrd):
        if verbose: print("Amber Parm files already exist. Skipping build phase to save time.")
        milestone.openmm.prmtop_filename = prmtop
        milestone.openmm.inpcrd_filename = inpcrd
        return # save us a little time by only running this when it matters

    leapfile = open(leapfilename,'w')
    leapfile.write(leap_string)
    leapfile.close()

    leapcmd = leap_program+' -f '+leapfilename+' > '+os.path.join(building,'leap.out')
    if verbose: print('running leap using following command:', leapcmd)

    errcode = os.system(leapcmd)
    # check to make sure everything is ok
    if errcode != 0: # then we hit a problem
        errormsg = "LEaP did not run properly. See %s/leap.out for details" % building
        raise Exception(errormsg)
    if (not os.path.exists(prmtop)) or (not os.path.exists(inpcrd)):
        errormsg = "LEaP did not generated expected prmtop & inpcrd files. See %s/leap.out for details" % building
        raise Exception(errormsg)
    milestone.openmm.prmtop_filename = prmtop
    milestone.openmm.inpcrd_filename = inpcrd
    return

def create_restraints(system, topology, constrained_list):
    '''Freeze receptor and ligand molecules, letting water and ions relax.
    Input:
     - system: the OpenMM system object for this milestone
     - topology: the OpenMM topology object of the MD simulation stages
     - constrained_list: a list of integers representing the indices of the atoms
         to constrain.
    Output:
     - None
    '''
    for i, atom in enumerate(topology.atoms()):
        if i in constrained_list:
            system.setParticleMass(i, 0*dalton)
    return


def create_simulation(seekrcalc, milestone):
    '''create the OpenMM simulation object.
    Input:
     - seekrcalc: The SeekrCalculation object that contains all the settings for
         the SEEKR calculation.
     - milestone: the Milestone() object to prepare the AMBER system for
    Output:
     - None
    '''
    if not milestone.openmm.prmtop_filename or not milestone.openmm.inpcrd_filename:
        print("Amber parameter/topology/coordinate files not found for milestone: %d. Skipping." % milestone.index)
        return
    prmtop_filename = milestone.openmm.prmtop_filename
    inpcrd_filename = milestone.openmm.inpcrd_filename
    prmtop = AmberPrmtopFile(prmtop_filename)
    inpcrd = AmberInpcrdFile(inpcrd_filename)

    system = prmtop.createSystem(nonbondedMethod=PME, nonbondedCutoff=1*nanometer) #, constraints=HBonds) # This is fine because h-bonds are always constrained in water!
    integrator = deepcopy(seekrcalc.min_equil.temp_equil_integrator) #LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)
    platform = seekrcalc.openmm.platform #Platform.getPlatformByName('CUDA')
    properties = seekrcalc.openmm.properties #{'CudaDeviceIndex':'0', 'CudaPrecision':'double'}

    # add restraints
    create_restraints(system, prmtop.topology, seekrcalc.min_equil.constrained) # TODO: clean this up
    barostat = MonteCarloBarostat(1.0*bar, seekrcalc.master_temperature*kelvin, 25)
    system.addForce(barostat)

    simulation = Simulation(prmtop.topology, system, integrator, platform, properties)
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    milestone.openmm.system = system
    milestone.openmm.simulation = simulation
    milestone.building_box_vectors = inpcrd.boxVectors
    return

def save_restart(seekrcalc, milestone, pdb_save_filename=None):
    '''Save an AMBER inpcrd file for easy restart.
    Input:
     - seekrcalc: The SeekrCalculation object that contains all the settings for
         the SEEKR calculation.
     - milestone: the Milestone() object to prepare the AMBER system for
    Output:
     - None
    '''
    if not milestone.openmm.prmtop_filename or not milestone.openmm.inpcrd_filename:
        print("Amber parameter/topology/coordinate files not found for milestone: %d. Skipping." % milestone.index)
        return
    state = milestone.openmm.simulation.context.getState(getPositions = True, enforcePeriodicBox = True)
    positions = state.getPositions()
    amber_parm = parmed.amber.AmberParm(milestone.openmm.prmtop_filename, milestone.openmm.inpcrd_filename)
    amber_parm.positions = positions
    if not pdb_save_filename:
        pdb_save_filename = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'temp_equil', 'equilibrated.pdb')
    amber_parm.save(pdb_save_filename, overwrite=True)
    return pdb_save_filename

