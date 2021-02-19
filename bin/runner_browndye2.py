"""
runner_browndye.py

A command-line tool for running SEEKR calculations using the BrownDye
engine either locally or on a computing resource.
"""
import os
import subprocess
import re
import argparse
import glob
import shutil
import xml.etree.ElementTree as ET
from xml.dom import minidom
from collections import defaultdict

import base
import seekr
#import openmmvt.base as base
import sim_browndye2
#import openmmvt.sim_browndye2 as sim_browndye2
REACTION_FILENAME = "rxns.xml"

def make_empty_pqrxml(directory, filename="empty.pqrxml"):
    """
    Create an empty pqrxml which is used in the Extraction phase to 
    generate new ligand encounter complex PQRs.
    """
    empty_pqrxml_path = os.path.join(directory, filename)
    with open(empty_pqrxml_path, "w") as f:
        empty_pqrxml_string = "<roottag> \n</roottag>"
        f.write(empty_pqrxml_string)
    return empty_pqrxml_path

def make_browndye_input_xml(model, rootdir, receptor_xml_filename, 
                            ligand_xml_filename, num_bd_steps, 
                            bd_directory=None, make_apbs_mode=True):
    """
    
    If bd_directory is None, then b_surface is assumed
    """    
    root = sim_browndye2.Root()
    if bd_directory is None:
        # B-surface
        bd_directory = os.path.join(
            rootdir, model.k_on_info.b_surface_directory)
        root.system.start_at_site = "False"
    else:
        # BD Milestone
        bd_directory = os.path.join(rootdir, bd_directory)
        root.system.start_at_site = "True"
        
    root.n_trajectories = num_bd_steps
    if root.n_trajectories_per_output > root.n_trajectories:
        root.n_trajectories_per_output = root.n_trajectories
    root.n_threads = model.browndye.num_threads
    root.system.reaction_file = "rxns.xml"
    reaction_filename = root.system.reaction_file
    if make_apbs_mode:
        root.system.solvent.debye_length = -1.0
        input_xml_filename = os.path.join(bd_directory, "apbs_input.xml")
    else:
        root.system.solvent.debye_length = model.browndye_settings.debye_length
        assert root.system.solvent.debye_length > 0.0, "The Debye length must "\
            "be set if make_apbs_mode=False"
        input_xml_filename = os.path.join(bd_directory, "input.xml")
        
    root.system.solvent.kT = model.master_temperature / 298.0
    root.system.time_step_tolerances.minimum_core_dt = 0.2
    root.system.time_step_tolerances.minimum_core_reaction_dt = 0.05
    
    root.system.solvent.ions = []
    for model_ion in model.browndye.apbs.ions:
        bd_ion = sim_browndye2.Ion()
        bd_ion.radius = model_ion.radius
        bd_ion.charge = model_ion.charge
        bd_ion.conc = model_ion.concentration
        root.system.solvent.ions.append(bd_ion)
        
    receptor_group = sim_browndye2.Group()
    receptor_group.name = sim_browndye2.BROWNDYE_RECEPTOR
    receptor_core = sim_browndye2.Core()
    receptor_core.name = sim_browndye2.BROWNDYE_RECEPTOR
    receptor_pqrxml_filename = os.path.basename(receptor_xml_filename)
    receptor_core.atoms = receptor_pqrxml_filename
    receptor_core.grid_spacing = model.browndye.apbs.inputgen.resolution
    receptor_group.core_list.append(receptor_core)
    ligand_group = sim_browndye2.Group()
    ligand_group.name = sim_browndye2.BROWNDYE_LIGAND
    ligand_core = sim_browndye2.Core()
    ligand_core.name = sim_browndye2.BROWNDYE_LIGAND
    ligand_pqrxml_filename = os.path.basename(ligand_xml_filename)
    ligand_core.atoms = ligand_pqrxml_filename
    ligand_core.grid_spacing = model.browndye.apbs.inputgen.resolution
    ligand_group.core_list.append(ligand_core)
    root.system.group_list.append(receptor_group)
    root.system.group_list.append(ligand_group)
    if make_apbs_mode:
        root.write(input_xml_filename)
        debye_length = sim_browndye2.make_and_run_apbs(
            root, input_xml_filename, 
            new_input_xml_base=sim_browndye2.BROWNDYE_INPUT_FILENAME)
    else:
        root.write(input_xml_filename)
        debye_length = model.browndye_settings.debye_length
                
    return debye_length, reaction_filename
    
def make_browndye_reaction_xml(model, abs_reaction_path, receptor_pqr_filename, 
                               ligand_pqr_filename, bd_milestone=None):
    """
    
    If bd_directory is None, then b_surface is assumed
    """
    rxnroot = sim_browndye2.Reaction_root()
    if bd_milestone is None:
        rxnroot.first_state = "b_surface"
        bd_milestone_name = "b_surface"
        bd_milestone2 = model.milestones[-1]
    else:
        bd_milestone2 = model.milestones[-2]
        rxnroot.first_state = str(bd_milestone.index)
        bd_milestone_name = bd_milestone.name
        
    
    ghost_indices_rec = []
    ghost_indices_lig = []
    print("adding ghost atom to file:", receptor_pqr_filename)
    ghost_index_rec = \
        sim_browndye2.add_ghost_atom_to_pqr_from_atoms_center_of_mass(
            receptor_pqr_filename, model.selections.site_com_indices)
    print("adding ghost atom to file:", ligand_pqr_filename)
    ghost_index_lig = \
        sim_browndye2.add_ghost_atom_to_pqr_from_atoms_center_of_mass(
            ligand_pqr_filename, model.selections.lig_com_indices)
    ghost_indices_rec.append(ghost_index_rec)
    ghost_indices_lig.append(ghost_index_lig)
        
    #model.browndye_settings.ghost_indices_rec = ghost_indices_rec
    #model.browndye_settings.ghost_indices_lig = ghost_indices_lig
    
    #ghost_indices_rec = model.browndye_settings.ghost_indices_rec
    #ghost_indices_lig = model.browndye_settings.ghost_indices_lig
    
    
    ghost_index_rec = ghost_indices_rec[0] # comes from the model?
    ghost_index_lig = ghost_indices_lig[0]
    rxn = sim_browndye2.Reaction()
    pair = sim_browndye2.Pair()
    
    rxn.name = str(bd_milestone2.index)
    rxn.state_after = str(bd_milestone2.index)
    pair.distance = bd_milestone2.radius
        
    rxn.state_before = rxnroot.first_state
    rxn.molecule0_group = sim_browndye2.BROWNDYE_RECEPTOR
    rxn.molecule0_core = sim_browndye2.BROWNDYE_RECEPTOR
    rxn.molecule1_group = sim_browndye2.BROWNDYE_LIGAND
    rxn.molecule1_core = sim_browndye2.BROWNDYE_LIGAND
    rxn.n_needed = 1
    pair.atom1_index = ghost_index_rec
    pair.atom2_index = ghost_index_lig
    
    rxn.pair_list.append(pair)
    rxnroot.reaction_list.append(rxn)
    
    rxnroot.write(abs_reaction_path)
    return

def run_bd_top(browndye_bin_dir, bd_directory, force_overwrite=False):
    """
    
    """
    curdir = os.getcwd()
    print("moving to directory:", bd_directory)
    os.chdir(bd_directory)
    simulation_filename = sim_browndye2.BROWNDYE_RECEPTOR + "_" \
        + sim_browndye2.BROWNDYE_LIGAND + "_simulation.xml"
    if os.path.exists(simulation_filename):
        if force_overwrite:
            print("force_overwrite set to True: existing files will be "\
                  "overwritten.")
            os.remove(simulation_filename)
        else:
            print("This anchor already has existing output files and the "\
                  "entered command would overwrite them. If you desire to "\
                  "overwrite the existing files, then use the "\
                  "--force_overwrite (-f) option, and all outputs will be "\
                  "deleted and replaced by a new run.")
            raise Exception("Cannot overwrite existing Browndye outputs.")
    bd_command = os.path.join(browndye_bin_dir, "bd_top")
    command = bd_command + " " + sim_browndye2.BROWNDYE_INPUT_FILENAME
    assert os.path.exists(sim_browndye2.BROWNDYE_INPUT_FILENAME), \
        "Necessary file doesn't exist: %s" % \
        sim_browndye2.BROWNDYE_INPUT_FILENAME
    print("running command:", command)
    os.system(command)
    assert os.path.exists(simulation_filename), "Problem occurred running "\
        "bd_top: simulation file %s was not generated." % simulation_filename
    os.chdir(curdir)
    return

def modify_variables(bd_milestone_directory, n_trajectories, n_threads=1, seed=0,
                     output_file=None):
    """
    
    """
    simulation_filename_base = sim_browndye2.BROWNDYE_RECEPTOR + "_" \
        + sim_browndye2.BROWNDYE_LIGAND + "_simulation.xml"
    simulation_filename = os.path.join(bd_milestone_directory, 
                                       simulation_filename_base)
    sim_file_old_lines = []
    sim_file_new_lines = []
    with open(simulation_filename, 'r') as f:
        for line in f.readlines():
            sim_file_old_lines.append(line)
    
    for line in sim_file_old_lines:
        new_line = line
        if n_trajectories is not None:
            new_line = re.sub(r"(?is)<n_trajectories>.+</n_trajectories>", 
                   "<n_trajectories> %d </n_trajectories>" % n_trajectories, 
                   new_line)

            new_line = re.sub(r"(?is)<n_trajectories_per_output>.+"\
                              "</n_trajectories_per_output>", 
                   "<n_trajectories_per_output> %d "\
                   "</n_trajectories_per_output>" % n_trajectories, 
                   new_line)
        
        if n_threads is not None:
            new_line = re.sub(r"(?is)<n_threads>.+</n_threads>", 
                   "<n_threads> %d </n_threads>" % n_threads, new_line)
            
        if seed is not None:
            new_line = re.sub(r"(?is)<seed>.+</seed>", 
                   "<seed> %d </seed>" % seed, new_line)
            
        if output_file is not None:
            new_line = re.sub(r"(?is)<output>.+</output>", 
                   "<output> %s </output>" % output_file, new_line)
            
        sim_file_new_lines.append(new_line)
    
    with open(simulation_filename, 'w') as f:
        for line in sim_file_new_lines:
            f.write(line)
            
    return

def run_nam_simulation(browndye_bin_dir, bd_directory, bd_mmvt_output_glob):
    """
    
    """
    curdir = os.getcwd()
    print("moving to directory:", bd_directory)
    os.chdir(bd_directory)
    bd_command = os.path.join(browndye_bin_dir, "nam_simulation")
    simulation_filename = sim_browndye2.BROWNDYE_RECEPTOR + "_" \
        + sim_browndye2.BROWNDYE_LIGAND + "_simulation.xml"
    command = bd_command + " " + simulation_filename
    print("running command:", command)
    os.system(command)
    results_glob = glob.glob(bd_mmvt_output_glob)
    assert len(results_glob) > 0, "Problem occurred running "\
        "nam_simulation: results file was not generated."
    os.chdir(curdir)
    return

def make_proc_file_last_frame(input_filename, output_filename, 
                              pqrxml_path_1, pqrxml_path_2):
    """
    Extract the last frame from a process_trajectories output XML
    and write a new XML containing only the last frame.
    """
    input_tree = ET.parse(input_filename)
    output_trajectory = ET.Element("trajectory")
    output_trajectory.text = "\n  "
    input_trajectory = input_tree.getroot()
    last_state = None
    for item in input_trajectory:
        if item.tag == "state":
            last_state = item
        elif item.tag == "atom_files":
            new_atom_files = ET.SubElement(output_trajectory, "atom_files")
            new_atom_files.text = " %s %s " % (pqrxml_path_1, pqrxml_path_2)
            new_atom_files.tail = "\n  "
        else:
            output_trajectory.append(item)
    assert last_state is not None
    output_trajectory.append(last_state)
    xmlstr = ET.tostring(output_trajectory).decode("utf-8")
    with open(output_filename, "w") as f:
        f.write(xmlstr)
        
    return

def extract_bd_surface(model, bd_milestone, max_b_surface_trajs_to_extract,
                       force_overwrite=False):
    """
    
    TODO: create a function to generate the correct trajectory file
    for PQR extraction.
    
    """
    pqr_lig_filenames = []
    pqr_rec_filenames = []
    b_surface_dir = model.browndye.b_surface_path
    b_surface_ligand_pqr = model.browndye.lig_dry_pqr_filename
    b_surface_ligand_pqrxml = os.path.splitext(
        b_surface_ligand_pqr)[0] + ".xml"
    b_surface_ligand_pqrxml_full_path = os.path.join(b_surface_dir, 
                                                     b_surface_ligand_pqrxml)
    assert os.path.exists(b_surface_ligand_pqrxml_full_path), "PQRXML file %s "\
        "not found for b-surface." % b_surface_ligand_pqrxml_full_path
    b_surface_receptor_pqr = model.browndye.rec_dry_pqr_filename
    b_surface_receptor_pqrxml = os.path.splitext(
        b_surface_receptor_pqr)[0] + ".xml"
    b_surface_receptor_pqrxml_full_path = os.path.join(
        b_surface_dir, b_surface_receptor_pqrxml)
    assert os.path.exists(b_surface_receptor_pqrxml_full_path), "PQRXML file "\
        "%s not found for b-surface." % b_surface_receptor_pqrxml_full_path
    bd_milestone_directory = os.path.join(model.project.rootdir, 
                                          bd_milestone.directory)
    extract_directory = os.path.join(bd_milestone_directory, 
                                     "extracted_from_b_surface")
    if not os.path.exists(extract_directory):
        os.mkdir(extract_directory)
    sitename = str(bd_milestone.index)
    empty_pqrxml_path = make_empty_pqrxml(extract_directory)
    process_trajectories = "process_trajectories"
    vtf_trajectory = "vtf_trajectory"
    
    quitting = False
    counter = 0
    for i in range(model.browndye.num_threads):
        if quitting: break
        print("extracting trajectories from traj number:", i)
        output_filename = os.path.join(extract_directory, 
                                      "rxn_output%d.txt" % i)
        if os.path.exists(output_filename):
            if force_overwrite:
                print("force_overwrite set to True: existing files will be "\
                      "overwritten.")
            else:
                print("This folder already has existing output files and the "\
                      "entered command would overwrite them. If you desire to "\
                      "overwrite the existing files, then use the "\
                      "--force_overwrite (-f) option, and all outputs will be "\
                      "deleted and replaced by a new run.")
                raise Exception("Cannot overwrite existing Browndye outputs.")
        traj_filename = os.path.join(b_surface_dir, "traj%d.xml" % i)
        trajindex_filename = os.path.join(b_surface_dir, "traj%d.index.xml" % i)
        assert os.path.exists(traj_filename), "trajectory output file %s not "\
            "found for b-surface. Are you sure you ran b-surface simulations?" \
            % traj_filename
        assert os.path.exists(trajindex_filename), "trajectory output file %s "\
            "not found for b-surface. Are you sure you ran b-surface "\
            "simulations?" % traj_filename
        command = "echo 'Browndye Trajectory number'; "\
            +process_trajectories+" -traj %s -index %s -srxn %s > %s" \
            % (traj_filename, trajindex_filename, sitename, 
               output_filename)
        print("running command:", command)
        std_out = subprocess.check_output(command, stderr=subprocess.STDOUT,
                                           shell=True)
        assert os.path.exists(output_filename) and \
                os.stat(output_filename).st_size > 0.0, "Problem running "\
            "process_trajectories: reaction list file not generated."
        number_list = []
        subtraj_list = []
        with open(output_filename, "r") as f:
            for line in f.readlines():
                if re.search("<number>",line):
                    number_list.append(int(line.strip().split()[1]))
                elif re.search("<subtrajectory>",line):
                    subtraj_list.append(int(line.strip().split()[1]))
                    
        number_list.sort()
        assert len(number_list) > 0, "No trajectories found in b_surface "\
            "simulations. Consider using larger outermost milestone or "\
            "simulating more b-surface trajectories."
        for j, rxn_number in enumerate(number_list):
            print("max_b_surface_trajs_to_extract:", max_b_surface_trajs_to_extract)
            if counter > max_b_surface_trajs_to_extract:
                quitting = True
                break
            rxn_subtraj = subtraj_list[j]
            proc_traj_basename = os.path.join(extract_directory,
                                              "proc_traj%d_%d" % (i, j))
            xml_traj_filename = proc_traj_basename + ".xml"
            command = process_trajectories+" -traj %s -index %s -n %d -sn %d "\
                "-nstride 1 > %s" % (traj_filename, trajindex_filename, 
                                     rxn_number, rxn_subtraj, xml_traj_filename)
            print("running command:", command)
            std_out = subprocess.check_output(command, stderr=subprocess.STDOUT,
                                           shell=True)
            assert os.path.exists(xml_traj_filename) and \
                os.stat(xml_traj_filename).st_size > 0.0, "Problem running "\
            "process_trajectories: trajectory XML file not generated."
            
            last_frame_name = proc_traj_basename + "_last.xml"
            make_proc_file_last_frame(xml_traj_filename, last_frame_name,
                                      empty_pqrxml_path, 
                                      b_surface_ligand_pqrxml_full_path)
            
            # write the last frame as a pqr file
            pqr_filename = os.path.join(extract_directory, 
                                        "lig%d_%d.pqr" % (i,j))
            pqr_lig_filenames.append(pqr_filename)
            command = vtf_trajectory+" -traj %s -pqr > %s"\
                % (last_frame_name, pqr_filename)
            print("running command:", command)
            std_out = subprocess.check_output(command, stderr=subprocess.STDOUT,
                                           shell=True)
            assert os.path.exists(pqr_filename) and \
                os.stat(pqr_filename).st_size > 0.0, "Problem running "\
                "vtf_trajectory: ligand%d_%d PQR file not generated." % (i, j)
                
            make_proc_file_last_frame(xml_traj_filename, last_frame_name,
                                  b_surface_receptor_pqrxml_full_path, 
                                  empty_pqrxml_path)
            
            pqr_rec_filename = os.path.join(extract_directory, 
                                            "receptor%d_%d.pqr" % (i,j))
            pqr_rec_filenames.append(pqr_rec_filename)
            command = vtf_trajectory+" -traj %s -pqr > "\
                "%s" % (last_frame_name, pqr_rec_filename)
            print("running command:", command)
            std_out = subprocess.check_output(
                command, stderr=subprocess.STDOUT, shell=True)
            assert os.path.exists(pqr_filename) and \
                os.stat(pqr_filename).st_size > 0.0, "Problem running "\
                "vtf_trajectory: receptor%d_%d PQR file not generated." % (i,j)
            
            os.remove(xml_traj_filename)
            counter += 1
    
    return pqr_lig_filenames, pqr_rec_filenames
    
def make_fhpd_directories(model, bd_milestone, lig_pqr_filenames, 
                          rec_pqr_filenames):
    """
    
    """
    b_surface_dir = model.browndye.b_surface_path
    bd_milestone_directory = os.path.join(model.project.rootdir, 
                                          bd_milestone.directory)
    extract_directory = os.path.join(bd_milestone_directory, 
                                     "extracted_from_b_surface")
    fhpd_directory = os.path.join(bd_milestone_directory, 
                                  "first_hitting_point_distribution")
    if not os.path.exists(fhpd_directory):
        os.mkdir(fhpd_directory)
    sitename = str(bd_milestone.index)
    #lig_pqr_glob = os.path.join(extract_directory, "lig*.pqr")
    #lig_pqr_filenames = glob.glob(lig_pqr_glob)
    #lig_pqr_orig_filenames = base.order_files_numerically(lig_pqr_filenames)
    #rec_pqr_orig_filename = os.path.join(extract_directory, "receptor.pqr")
    #rec_pqr_orig_filenames = base.order_files_numerically(rec_pqr_filenames)
    
    directories = []
    for i, (lig_pqr_orig_filename, rec_pqr_orig_filename) in enumerate(zip(
            lig_pqr_filenames, rec_pqr_filenames)):
        ligand_basename = os.path.splitext(
            os.path.basename(lig_pqr_orig_filename))[0]
        directory_name = os.path.join(fhpd_directory, ligand_basename)
        if not os.path.exists(directory_name):
            os.mkdir(directory_name)
        
        receptor_xml_filename = os.path.join(directory_name, 'receptor.xml')
        rec_pqr_filename = os.path.join(directory_name, "receptor.pqr")
        shutil.copyfile(rec_pqr_orig_filename, rec_pqr_filename)
        ligand_xml_filename = os.path.join(directory_name, 'ligand.xml')
        lig_pqr_filename = os.path.join(directory_name, "ligand.pqr")
        shutil.copyfile(lig_pqr_orig_filename, lig_pqr_filename)
        
        debye_length, reaction_filename = make_browndye_input_xml(
            model, model.project.rootdir, receptor_xml_filename, 
            ligand_xml_filename, model.browndye.fhpd_numtraj, 
            bd_directory=directory_name, make_apbs_mode=True)
        abs_reaction_path = os.path.join(directory_name, 
                                         reaction_filename)
        make_browndye_reaction_xml(model, abs_reaction_path, rec_pqr_filename, 
                                   lig_pqr_filename, bd_milestone=None)
        
        sim_browndye2.make_pqrxml(rec_pqr_filename, 
            browndye2_bin="", output_xml_filename=receptor_xml_filename)
        
        sim_browndye2.make_pqrxml(lig_pqr_filename, 
            browndye2_bin="", output_xml_filename=ligand_xml_filename)
        local_receptor_xml_filename = os.path.join(
            directory_name, 'receptor.xml')
        #shutil.copyfile(receptor_xml_filename, local_receptor_xml_filename)
        directories.append(directory_name)
        
    return directories

def combine_fhpd_results(model, bd_milestone, fhpd_directories):
    """
    
    """
    reaction_dict = defaultdict(float)
    number_escaped = 0
    number_stuck = 0
    number_total = 0
    number_total_check = 0
    results_filename_list = []
    for fhpd_directory in fhpd_directories:
        results_glob = os.path.join(fhpd_directory, 
                                    "results.xml")
        results_filename_list += glob.glob(results_glob)
    
    assert len(results_filename_list) > 0, "No BD output files found."
    for results_filename in results_filename_list:
        tree = ET.parse(results_filename)
        root = tree.getroot()
        reactions_XML = root.find("reactions")
        number_total += int(reactions_XML.find("n_trajectories").text.strip())
        number_stuck += int(reactions_XML.find("stuck").text.strip())
        number_escaped += int(reactions_XML.find("escaped").text.strip())
        for completed_XML in reactions_XML.iter("completed"):
            name = completed_XML.find("name").text.strip()
            n = int(completed_XML.find("n").text.strip())
            reaction_dict[name] += n
            number_total_check += n
        
    assert number_total == number_total_check + number_stuck + number_escaped
    for completed_XML in reactions_XML.iter("completed"):
        reactions_XML.remove(completed_XML)
    
    reactions_XML.find("n_trajectories").text = str(number_total)
    reactions_XML.find("stuck").text = str(number_stuck)
    reactions_XML.find("escaped").text = str(number_escaped)
    for key in reaction_dict:
        completed_XML = ET.SubElement(reactions_XML, "completed")
        completed_XML.text = "\n      "
        completed_XML.tail = "\n  "
        name_XML = ET.SubElement(completed_XML, "name")
        name_XML.text = key
        name_XML.tail = "\n      "
        n_XML = ET.SubElement(completed_XML, "n")
        n_XML.text = str(int(reaction_dict[key]))
        n_XML.tail = "\n    "
    
    xmlstr = ET.tostring(root).decode("utf-8")
    bd_milestone_directory = os.path.join(model.project.rootdir, 
                                          bd_milestone.directory)
    dest_filename = os.path.join(bd_milestone_directory, "results.xml")
    with open(dest_filename, 'w') as f:
        f.write(xmlstr)
        
    return

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "bd_milestone", metavar="BD_MILESTONE", type=str, 
        help="Which milestone to run BrownDye for. Arguments may be the "\
        "string 'b_surface' or a numerical index of a given BD_milestone.")
    argparser.add_argument(
        "input_file", metavar="INPUT_FILE", type=str, 
        help="name of input file for OpenMMVT calculation. This would be the "\
        "XML file generated in the prepare stage.")
    """ # TODO: remove
    argparser.add_argument("-r", "--run", dest="run", default=False,
        help="Run the BrownDye bd_top and nam_simulation programs for the "\
        "BD_MILESTONE indicated. By default, this mode is activated if the "\
        "extract argument is False, otherwise, if extract is True, this "\
        "argument must be specified in order to do both.",
        action="store_true")
    argparser.add_argument("-e", "--extract", dest="extract", default=False,
        help="Whether to perform an extraction of starting positions from "\
        "b-surface simulations. This mode is incompatible with an argument of "\
        "'b_surface' for argument BD_MILESTONE", action="store_true")
    """
    argparser.add_argument("-n", "--n_trajectories", dest="n_trajectories", 
        default=None, help="Enter a different number of trajectories to run "\
        "the simulation if a different number of trajectories are desired "\
        "than what is in the input.xml file.", type=int)
    argparser.add_argument("-t", "--n_threads", dest="n_threads", 
        default=None, help="Enter a different number of threads to run "\
        "the simulation if a different number of threads are desired "\
        "than what is in the input.xml file.", type=int)
    argparser.add_argument("-o", "--output_file", dest="output_file",
        default=None, help="Enter a new output file name different from the "\
        "output file specified in the input.xml file.", type=str)
    argparser.add_argument("-s", "--seed", dest="seed", default=False, 
        help="Enter a new random number seed if different from the "\
        "seed specified in the input.xml file.", type=int)
    argparser.add_argument("-m", "--max_b_surface_trajs", 
        dest="max_b_surface_trajs", help="In extraction mode, enter the "\
        "maximum number of successful encounter states to extract from the "\
        "b-surface simulations to start the bd_milestone simulations from.", 
        type=int, default=1000)
    argparser.add_argument("-f", "--force_overwrite", dest="force_overwrite",
        default=False, help="Toggle whether to overwrite existing files "\
        "within milestone. If this option is enabled, "\
        "then the anchor's production directory will be emptied of all "\
        "output, trajectory, and results files for the new simulation.", 
        action="store_true")
        
    args = argparser.parse_args() # parse the args into a dictionary
    args = vars(args)
    bd_milestone_index = args["bd_milestone"]
    input_file = args["input_file"]
    #run = args["run"]
    #extract = args["extract"]
    n_trajectories = args["n_trajectories"]
    n_threads = args["n_threads"]
    output_file = args["output_file"]
    seed = args["seed"]
    max_b_surface_trajs_to_extract = args["max_b_surface_trajs"]
    force_overwrite = args["force_overwrite"]
    
    b_surface = False
    if bd_milestone_index.lower() == "b_surface":
        bd_milestone_index = bd_milestone_index.lower()
        b_surface = True
    else:
        try:
            bd_milestone_index = int(bd_milestone_index)
        except ValueError:
            print("The bd_milestone argument must either be the string "\
                  "'b_surface' or an integer.")
            exit()
    
    assert os.path.exists(input_file), "A nonexistent input file was provided."
    
    me = seekr.openSeekrCalc(input_file)
    
    if n_threads is not None:
        me.browndye.num_threads = n_threads
        
    assert os.path.exists(me.project.rootdir), "An incorrect directory "\
        "was provided."
    
    if b_surface:
        bd_milestone_directory = me.browndye.b_surface_path
        
    else:
        assert bd_milestone_index >= 0, "only positive indices allowed."
        bd_milestone = me.milestones[bd_milestone_index]
        assert bd_milestone.bd, "BD for this milestone not enabled."
        
        bd_milestone_directory = os.path.join(
            me.project.rootdir, bd_milestone.directory)
    
    if not b_surface:
        #assert not b_surface, "Extraction may not be performed on the "\
        #    "b-surface."
        pqr_lig_filenames, pqr_rec_filenames = extract_bd_surface(
            me, bd_milestone, max_b_surface_trajs_to_extract, 
            force_overwrite)
        bd_directory_list = make_fhpd_directories(
            me, bd_milestone, pqr_lig_filenames, pqr_rec_filenames)
        
    else:
        #run = True
        bd_directory_list = [bd_milestone_directory]
        #receptor_xml_filename = os.path.join(bd_milestone_directory, 
        #                                     'receptor.xml')
        #ligand_xml_filename = os.path.join(bd_milestone_directory, 'ligand.xml')
        rec_filename = os.path.join(bd_milestone_directory, "receptor.pqr")
        lig_filename = os.path.join(bd_milestone_directory, "ligand.pqr")
        receptor_xml_filename = os.path.join(bd_milestone_directory, "receptor.xml")
        ligand_xml_filename = os.path.join(bd_milestone_directory, "ligand.xml")
        shutil.copyfile(me.browndye.rec_dry_pqr_filename, rec_filename)
        shutil.copyfile(me.browndye.lig_dry_pqr_filename, lig_filename)
        
        abs_reaction_path = os.path.join(bd_milestone_directory, 
                                         REACTION_FILENAME)
        make_browndye_reaction_xml(me, abs_reaction_path, rec_filename, 
                                   lig_filename, bd_milestone=None, )
        receptor_xml_filename_dummy = sim_browndye2.make_pqrxml(rec_filename)
        ligand_xml_filename_dummy = sim_browndye2.make_pqrxml(lig_filename)
        
        debye_length, reaction_filename = make_browndye_input_xml(
            me, me.project.rootdir, receptor_xml_filename, 
            ligand_xml_filename, me.browndye.prods_per_anchor, 
            bd_directory=None, make_apbs_mode=True)
        
    for bd_directory in bd_directory_list:
        print("bd_directory:", bd_directory)
        run_bd_top("", bd_directory, force_overwrite)
        modify_variables(bd_directory, n_trajectories, n_threads,
                         seed, output_file)
        run_nam_simulation("", bd_directory, "results.xml")
    
    if not b_surface:
        combine_fhpd_results(me, bd_milestone, bd_directory_list)