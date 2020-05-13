'''
Created on June 7, 2018

@author: lvotapka

Functions and objects for running Umbrella sampling
'''

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
import os, time, glob, re
from sys import stdout
from seekr import amber
from base import saveStateWithoutParam
import mdtraj

verbose = True

def load_last_mdtraj_frame(dcd_filename, prmtop_filename, atom_indices=None):
    '''
    This function returns the last frame of a DCD file as an MDtraj object.
    Input:
     - dcd_filename: string that represents the path to the dcd file to load
     - prmtop_filename: string that represents the path to the parm7 file to load
    Output:
     - lastframe: an mdtraj Trajectory object that contains a single frame: the
       last one in the dcd.
    '''
    if atom_indices is not None:
        mytraj_iter = mdtraj.iterload(dcd_filename, top=prmtop_filename, atom_indices=atom_indices)
    else:
        mytraj_iter = mdtraj.iterload(dcd_filename, top=prmtop_filename) # Trajectory object

    for frame in mytraj_iter:
        lastframe = frame[-1]
    return lastframe

def create_forces(seekrcalc, milestone, system):
    '''
    Add the umbrella force: which maintains the ligand on the surface of the
    spherical milestone.
    Input:
     - seekrcalc: The SeekrCalculation object that contains all the settings for
         the SEEKR calculation.
     - milestone: the Milestone() object to run the simulation for
     - system: the OpenMM system object to add the force to
    Output:
     - None
    '''
    #new_force = CustomCentroidBondForce(2, '0.5*k*(z2-z1-length)^2')
    new_force = CustomCentroidBondForce(2, '0.5*k*(distance(g1,g2)-radius)^2')
    k = new_force.addGlobalParameter(
        'k', seekrcalc.umbrella_stage.force_constant)
    #r0 = new_force.addGlobalParameter('length', milestone.radius*angstrom)
    r0 = new_force.addGlobalParameter('radius', milestone.radius*angstrom)
    assert len(milestone.atom_selection_1) > 0
    g1 = new_force.addGroup(milestone.atom_selection_1)
    assert len(milestone.atom_selection_2) > 0
    g2 = new_force.addGroup(milestone.atom_selection_2)
    if verbose: print(("k:", seekrcalc.umbrella_stage.force_constant, 
                       "radius:", milestone.radius*angstrom, 
                       "g1:", milestone.atom_selection_1, 
                       "g2:", milestone.atom_selection_2))
    new_force.addBond([g1, g2], [])
    if verbose: print(("new_force.getNumGlobalParameters():", 
                       new_force.getNumGlobalParameters()))
    if verbose: print(("new_force.getNumPerBondParameters():", 
                       new_force.getNumPerBondParameters()))
    return new_force

def launch_umbrella_stage(seekrcalc, milestone, box_vectors=None, traj_name='umbrella1.dcd'):
    '''launch an umbrella sampling job.
    Input:
     - seekrcalc: The SeekrCalculation object that contains all the settings for
         the SEEKR calculation.
     - milestone: the Milestone() object to run the simulation for
    Output:
     - ending_box_vectors: An array of three vectors that represents the ending
     state of the periodic box. This may not be the same as was started, but could
     have changed through the course of a constant pressure simulation.
    '''
    prmtop_filename = milestone.openmm.prmtop_filename
    pdb_filename = milestone.openmm.umbrella_pdb_filename
    root_name = os.path.splitext(pdb_filename)[0]
    prev_end_state_filename = os.path.join(
        seekrcalc.project.rootdir, milestone.directory, 'md', 'umbrella', 
        '%s_end.state' % root_name)
    
    inpcrd_filename = milestone.openmm.inpcrd_filename
    if verbose: print("opening files:", prmtop_filename, inpcrd_filename, pdb_filename)
    prmtop = AmberPrmtopFile(prmtop_filename)
    inpcrd = AmberInpcrdFile(inpcrd_filename)
    
    # This is fine because h-bonds are always constrained in water!
    system = prmtop.createSystem(nonbondedMethod=PME, 
                                 nonbondedCutoff=1*nanometer, 
                                 constraints=HBonds) 
    integrator = LangevinIntegrator(seekrcalc.master_temperature*kelvin, 
                                    1/picosecond, 0.002*picoseconds)
    platform = Platform.getPlatformByName('CUDA')

    # TODO: change this back
    properties = seekrcalc.openmm.properties #{'CudaDeviceIndex':'0', 'CudaPrecision':'mixed', 'UseCpuPme':'false'}
    
    if os.path.exists(prev_end_state_filename):
        print('loaded state from', prev_end_state_filename)
    elif os.path.exists(pdb_filename):
        pdb = PDBFile(pdb_filename)
        my_positions = pdb.positions
        print('loading positions from', pdb_filename)
    else:
        basename = os.path.basename(pdb_filename)
        no_ext = os.path.splitext(basename)[0]
        dcd_filename = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'umbrella', '%s.dcd' % no_ext)
        assert os.path.exists(dcd_filename), "Cannot load DCD or PDB file for umbrella stage, none exist:" + dcd_filename
        print(("Restarting failed umbrella stage from last frame of DCD file: ", dcd_filename))
        last_fwd_frame = load_last_mdtraj_frame(dcd_filename, milestone.openmm.prmtop_filename)
        my_positions = last_fwd_frame.xyz[0]
        print('loading positions from', dcd_filename)

        # add restraints
    new_force = create_forces(seekrcalc, milestone, system) #system, prmtop.topology, inpcrd.positions, seekrcalc.min_equil.constrained)
    system.addForce(new_force)
    barostat_force_id = None
    if seekrcalc.umbrella_stage.barostat:
        barostat = MonteCarloBarostat(seekrcalc.umbrella_stage.barostat_pressure, seekrcalc.master_temperature*kelvin, seekrcalc.umbrella_stage.barostat_freq)
        barostat_force_id = system.addForce(barostat)
    
    simulation = Simulation(prmtop.topology, system, integrator, 
                                     platform, properties)
    
    if os.path.exists(prev_end_state_filename):
        simulation.loadState(prev_end_state_filename)
    else:
        simulation.context.setPositions(my_positions)
        simulation.context.setVelocitiesToTemperature(seekrcalc.master_temperature*kelvin)
        if box_vectors:
            simulation.context.setPeriodicBoxVectors(*box_vectors)
            print('box_vectors assigned', box_vectors)
        elif milestone.umbrella_box_vectors:
            simulation.context.setPeriodicBoxVectors(milestone.umbrella_box_vectors)
            print('milestone.umbrella_box_vectors assigned', box_vectors)
        elif inpcrd.boxVectors is not None:
            simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
            print('inpcrd.boxVectors assigned', box_vectors)
    
        if verbose: print(("Running energy minimization on milestone:", milestone.index))
        simulation.minimizeEnergy()

    umbrella_traj = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'umbrella', traj_name)
    simulation.reporters.append(StateDataReporter(stdout, seekrcalc.umbrella_stage.energy_freq, step=True, potentialEnergy=True, temperature=True, volume=True))
    simulation.reporters.append(DCDReporter(umbrella_traj, seekrcalc.umbrella_stage.traj_freq))
    starttime = time.time()

    state_filename = os.path.join(seekrcalc.project.rootdir, milestone.directory, 'md', 'umbrella', 'backup.state')
    step_chunk_size = min(seekrcalc.umbrella_stage.steps, 
                          seekrcalc.umbrella_stage.traj_freq)
    current_step = 0
    simulation.currentStep = 0
    print(("running %d steps" % seekrcalc.umbrella_stage.steps))
    while current_step < seekrcalc.umbrella_stage.steps:
        try:
            simulation.saveState(state_filename)
            simulation.step(step_chunk_size)
            current_step = current_step + step_chunk_size
        except Exception:
            print("Alert! NaN error detected. Restarting from saved state.")
            simulation.loadState(state_filename)
            simulation.currentStep = current_step

    #simulation.step(seekrcalc.umbrella_stage.steps) # old way: simulate steps directly
    print(("time:", time.time() - starttime, "s"))
    if barostat_force_id is not None:
        system.removeForce(barostat_force_id)
    end_state = simulation.context.getState(getPositions=True, 
                                            getVelocities=True)
    root_name = os.path.splitext(umbrella_traj)[0]
    end_state_filename = os.path.join(
        seekrcalc.project.rootdir, milestone.directory, 'md', 'umbrella', 
        '%s_end.state' % root_name)
    #simulation.saveState(end_state_filename)
    saveStateWithoutParam(my_simulation=simulation, file=end_state_filename)
    ending_box_vectors = end_state.getPeriodicBoxVectors()
    milestone.openmm.simulation = simulation
    milestone.umbrella_box_vectors = ending_box_vectors
    return ending_box_vectors, umbrella_traj

def generate_umbrella_filenames(seekr_calc, milestone):
    umbrella_file_glob = os.path.join(seekr_calc.project.rootdir, milestone.directory, 'md', 'umbrella', 'umbrella*.dcd')
    existing_umbrella_files = glob.glob(umbrella_file_glob)
    if not existing_umbrella_files: # then the directory is empty, we are starting over
        milestone.openmm.umbrella_pdb_filename = os.path.join(seekr_calc.project.rootdir, milestone.directory, 'md', 'temp_equil', 'equilibrated.pdb')
        new_dcd_filename = 'umbrella1.dcd'
        new_pdb_filename = 'umbrella1.pdb'
        milestone.umbrella_box_vectors = None
    else: # then some already exist
        number_list = []
        for existing_file in existing_umbrella_files:
            number_list.append(int(re.findall(r".+(\d+).dcd", existing_file)[0]))
        current_num = max(number_list)
        milestone.openmm.umbrella_pdb_filename = os.path.join(seekr_calc.project.rootdir, milestone.directory, 'md', 'umbrella', 'umbrella%d.pdb' % current_num)
        next_num = current_num + 1
        new_dcd_filename = 'umbrella%d.dcd' % next_num
        new_pdb_filename = 'umbrella%d.pdb' % next_num

    return new_dcd_filename, new_pdb_filename
