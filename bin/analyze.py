'''
Created on May 12, 2020

As OpenSEEKR performs its calculations, this module provides a number of tools
to analyze its transitions and trajectories.

Most notably, this module computes the rate constants and thermodynamics of
binding and unbinding for the ligand-receptor system.

@author: lvotapka
'''
import os
import sys
import glob
import xml.etree.ElementTree as ET
from collections import defaultdict
import random
import subprocess

import numpy as np
from numpy import log, exp
import scipy.linalg as la
import mdtraj
from scipy.stats import gamma
import matplotlib.pyplot as plt

import seekr
import dig_deeper

MATRIX_EXPONENTIAL = 9999999
GAS_CONSTANT = 0.0019872 # in kcal/mol*K
BOLTZMAN = 1.3806488e-23
DEFAULT_IMAGE_DIR = "images_and_plots/"


def quadriture(err1, err2):
    """
    Add two errors in quadriture.
    """
    return float(np.sqrt(err1**2 + err2**2))

def pretty_string_value_error(value, error, error_digits=2, use_unicode=True):
    """
    Returns a value/error combination of numbers in a scientifically 
    'pretty' format.
    
    Scientific quantities often come as a *value* (the actual 
    quantity) and the *error* (the uncertainty in the value).
        
    Given two floats, value and error, return the two in a 
    'pretty' formatted string: where the value and error are truncated
    at the correct precision.
    
    Parameters
    ----------
    value : float
        The quantity in question
        
    error : float
        The uncertainty of the quantity
        
    error_digits : int, default 2
        How many significant figures the error has. Scientific 
        convention holds that errors have 1 or (at most) 2 significant
        figures. The larger number of digits is chosen here by default.
        
    Returns
    -------
    new_string : str
        A new list of strings sorted numerically
    
    Examples
    --------
    
        >>> pretty_string_value_error(1.23456789e8, 4.5678e5, 
                                      error_digits=2)
        "1.2346 +/- 0.0046 * 10^+08"
    
        >>> pretty_string_value_error(5.6e-2, 2.0e-3, error_digits=1)
        "5.6 +/- 0.2 * 10^-02"
    
    
    """
    if error is None:
        if use_unicode:
            new_string = "{:.6E} \u00B1 UNKNOWN ERROR".format(value)
        else:
            new_string = "{:.6E} +/- UNKNOWN ERROR".format(value)
    else:
        if not np.isfinite(value):
            return str(value)
        assert "e" in "{:e}".format(value), "Cannot convert into scientific "\
            "notation: {1}".format(value)
        value_mantissa_str, value_exponent_str = \
            "{:e}".format(value).strip().split('e')
        value_mantissa = float(value_mantissa_str)
        value_exponent = int(value_exponent_str)
        error_mantissa_str, error_exponent_str = \
            "{:e}".format(error).strip().split('e')
        error_mantissa = float(error_mantissa_str)
        error_exponent = int(error_exponent_str)
        padding = value_exponent - error_exponent + error_digits - 1
        if padding < 1: padding = 1
        exp_diff = error_exponent - value_exponent
        string_for_formatting = "{:.%df}" % padding
        new_value_mantissa = string_for_formatting.format(value_mantissa)
        new_error_mantissa = string_for_formatting.format(
            error_mantissa*10**exp_diff)
        
        if use_unicode:
            new_string = "%s \u00B1 %s * 10^%s" % (
                new_value_mantissa, new_error_mantissa, value_exponent_str)
        else:
            new_string = "%s +/- %s * 10^%s" % (
                new_value_mantissa, new_error_mantissa, value_exponent_str)
    return new_string

def get_umbrella_avg_distance(me, milestone):
    umbrella_dir = os.path.join(me.project.rootdir, milestone.directory, 
                                'md', 'umbrella')
    prmtop = os.path.join(me.project.rootdir, milestone.directory, 'md', 
                          'building', 'holo.parm7')
    umbrella_dcd_list = glob.glob(os.path.join(umbrella_dir, 'umbrella*.dcd'))
    
    mytraj1 = mdtraj.load(umbrella_dcd_list, top=prmtop, 
                           atom_indices=milestone.atom_selection_1)
    mytraj2 = mdtraj.load(umbrella_dcd_list, top=prmtop, 
                           atom_indices=milestone.atom_selection_2)
    com1_array = mdtraj.compute_center_of_mass(mytraj1)
    com2_array = mdtraj.compute_center_of_mass(mytraj2)
    distance_total = 0.0
    for com1, com2 in zip(com1_array, com2_array):
        this_distance = np.linalg.norm(com2-com1)
        distance_total += this_distance
    
    avg_distance = 10.0 * distance_total / com1_array.shape[0]
    return avg_distance

def get_fwd_rev_avg_distance(me, milestone):
    fwd_rev_dir = os.path.join(me.project.rootdir, milestone.directory, 
                                'md', 'fwd_rev')
    prmtop = os.path.join(me.project.rootdir, milestone.directory, 'md', 
                          'building', 'holo.parm7')
    fwd_dcd_list = sorted(glob.glob(
        os.path.join(fwd_rev_dir, 'forward*.dcd')), 
        key=seekr.sort_forward_dcd_key)
    data_file_name = os.path.join(fwd_rev_dir, 'transition_fwd.dat')
    
    downward_indices = dig_deeper.read_data_file_transitions_down(
        data_file_name, destination='1')
    upward_indices = dig_deeper.read_data_file_transitions_down(
        data_file_name, destination='3')
    
    upward_distance_total = 0.0
    upward_counter = 0
    downward_distance_total = 0.0
    downward_counter = 0
    for i, dcd_file in enumerate(fwd_dcd_list):
        last_frame1 = seekr.load_last_mdtraj_frame(
            dcd_file, prmtop, atom_indices=milestone.atom_selection_1)
        last_frame2 = seekr.load_last_mdtraj_frame(
            dcd_file, prmtop, atom_indices=milestone.atom_selection_2)
        com1_array = mdtraj.compute_center_of_mass(last_frame1)
        com2_array = mdtraj.compute_center_of_mass(last_frame2)
        #print('com1_array:', com1_array, 'com2_array:', com2_array)
        this_distance = np.linalg.norm(com2_array[0]-com1_array[0])
        
        if i in upward_indices:
            upward_distance_total += this_distance
            upward_counter += 1
        elif i in downward_indices:
            downward_distance_total += this_distance
            downward_counter += 1
        if upward_counter >= 10 and downward_counter >= 10:
            # enough statistics
            break
    if upward_counter > 0:
        upward_avg_distance = 10.0 * upward_distance_total / upward_counter
    else:
        upward_avg_distance = None
        
    if downward_counter > 0:
        downward_avg_distance = 10.0 * downward_distance_total / \
            downward_counter
    else:
        downward_avg_distance = None
        
    return upward_avg_distance, downward_avg_distance

def deserealize_transition_info(seekrcalc, milestone):
    """
    
    """
    N_i_j_dict = defaultdict(int)
    transition_info_pickle = os.path.join(
        seekrcalc.project.rootdir, milestone.directory, 'md', 'fwd_rev', 
        'transition_info.xml')
    assert os.path.exists(transition_info_pickle), \
        "statistics not found for milestone: %d" % milestone.index
    tree = ET.parse(transition_info_pickle)
    root = tree.getroot()
    
    milestone_index = root.find('milestone_index').text
    milestone_siteid = root.find('milestone_siteid').text
    milestone_radius = root.find('milestone_radius').text
    avg_incubation_time = float(root.find('avg_incubation_time').text)
    for transition in root.find('transitions'):
        source = int(transition.find("source").text)
        destination = int(transition.find("destination").text)
        destRadius = transition.find("destRadius").text
        count = int(transition.find("count").text)
        N_i_j_dict[(source, destination)] += count
        
    return N_i_j_dict, avg_incubation_time

def browndye_run_compute_rate_constant(compute_rate_constant_program,
                                       results_filename, 
                                       sample_error_from_normal=False):
    """
    run the BrownDye program compute_rate_constant to find k_ons
    and the value k(b).
    
    Parameters:
    -----------
    compute_rate_constant_program : str
        The exact command to use to run the Browndye 
        compute_rate_constant program.
    
    results_filename : str
        A path to the results XML file, which is an output from
        running the BrownDye program.
    
    sample_error_from_normal : bool, default False
        Add a fluctuation to the k-on value(s) sampled from a
        normal distribution with a standard deviation equal to the
        k-on's error. This is used by the Monte Carlo error
        estimator to incorporate the error in the k_b value
        obtained from BrownDye.
    
    Results:
    --------
    k_ons : dict
        dictionary whose keys are various milestoning states and
        whose values are the k-ons to those states.
    
    k_on_errors : dict
        dictionary whose keys are various milestoning states and
        whose values are the errors in the k-ons to those states.
    
    reaction_probabilities : dict
        dictionary whose keys are various milestoning states and
        whose values are the probabilities of reaching those states
        from the b-surface.
    
    reaction_probability_errors : dict
        dictionary whose keys are various milestoning states and
        whose values are the errors in the probabilities of 
        reaching those states from the b-surface.
    
    """
    cmd = "%s < %s" % (compute_rate_constant_program, 
                       results_filename)
    #print("running command:", cmd)
    output_string = subprocess.check_output(cmd, shell=True)
    root = ET.fromstring(output_string)
    k_ons = {}
    k_on_errors = {}
    reaction_probabilities = {}
    reaction_probability_errors = {}
    k_b = None
    for rate in root.iter("rate"):
        name = int(rate.find("name").text.strip())
        rate_constant = rate.find("rate_constant")
        k_on_tag = rate_constant.find("mean")
        k_on_tag_high = rate_constant.find("high")
        assert k_on_tag is not None
        assert k_on_tag_high is not None
        k_on = float(k_on_tag.text)
        k_on_error = 0.5*(float(k_on_tag_high.text) - k_on)
        reaction_probability = rate.find("reaction_probability")
        beta_tag = reaction_probability.find("mean")
        beta_tag_high = reaction_probability.find("high")
        assert beta_tag is not None
        assert beta_tag_high is not None
        beta = float(beta_tag.text)
        beta_error = float(beta_tag_high.text) - beta
        reaction_probabilities[name] = beta
        reaction_probability_errors[name] = beta_error
        if beta != 0.0:
            k_b = k_on / beta
        
        if sample_error_from_normal:
            k_on = np.random.normal(loc=k_on, scale=k_on_error)
            
        k_ons[name] = k_on
        k_on_errors[name] = k_on_error
    
    assert k_b is not None, "No BD reactions from the b-surface " \
        "successfully reached any of the milestone surfaces."
    #k_ons["b_surface"] = k_b
    #k_on_errors["b_surface"] = 0.0
    #reaction_probabilities["b_surface"] = 1.0
    #reaction_probability_errors["b_surface"] = 0.0
    return k_on, k_on_error, beta, beta_error
        
def browndye_parse_bd_milestone_results(results_filename, 
                                        sample_error_from_normal=False):
    """
    Read and extract transition probabilities for a BD milestone.
    
    Parameters:
    -----------
    results_filename : str
        A path to the results XML file for a BD milestone, which 
        is an output from running the BrownDye program.
        
    sample_error_from_normal : bool, default False
        Add a fluctuation to the probabilities sampled from a
        normal distribution with a standard deviation equal to
        1/sqrt(n-1). This is used by the Monte Carlo error
        estimator to incorporate the error in the probabilities.
    
    Results:
    --------
    transition_probabilities : dict
        The probabilities of transitions.
    """
    transition_counts = {}
    transition_probabilities = {}
    assert os.path.exists(results_filename), "You must perform successful "\
        "extractions and runs of all BD milestones if k-on settings are "\
        "provided in the model XML. Missing file: " \
        + results_filename
    root = ET.parse(results_filename)
    reactions = root.find("reactions")
    n_trajectories = int(reactions.find("n_trajectories").text.strip())
    stuck = int(reactions.find("stuck").text.strip())
    escaped = int(reactions.find("escaped").text.strip())
    completed_prob = 0.0
    completed_count = 0
    for completed in root.iter("completed"):
        name = int(completed.find("name").text.strip())
        n = int(completed.find("n").text.strip())
        transition_counts[name] = n
        completed_count += n
        
    transition_counts["escaped"] = escaped
    for key in transition_counts:
        n = transition_counts[key]
        avg = n / (completed_count + escaped)
        if sample_error_from_normal:
            assert n > 1, "Too few transitions to compute error."
            std_dev = avg / np.sqrt(n-1)
            transition_probabilities[key] = np.random.normal(
                loc=avg, scale=std_dev)
            
        else:
            transition_probabilities[key] = avg
            
        completed_prob += transition_probabilities[key]
        
    transition_probabilities["escaped"] = 1.0 - completed_prob
    return transition_probabilities

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

def make_matrices(me):
    """
    
    """
    n = len(me.milestones)
    K = np.zeros((n,n), dtype=np.double)
    N_ij = np.zeros((n,n), dtype=np.double)
    K_hat = np.zeros((n-1,n-1), dtype=np.double)
    avg_t = np.zeros((n-1,1), dtype=np.double)
    R_i = np.zeros((n,1), dtype=np.double)
    bd_milestone_prob = 0.0
    
    for milestone in me.milestones:
        if milestone.md:
            N_i_j_dict, avg_inc_time = deserealize_transition_info(me, milestone)
            avg_t[milestone.index, 0] = avg_inc_time
            N_i = 0.0
            for key in N_i_j_dict:
                source = key[0]
                dest = key[1]
                N_i += N_i_j_dict[key]
                
            for key in N_i_j_dict:
                source = key[0]
                dest = key[1]
                if source == 0:
                    K[0,0] = 1.0
                    K_hat[0,1] = 1.0
                    N_ij[0,dest] = N_i_j_dict[key]
                    R_i[0] = avg_t[source,0]
                else:
                    K[source, dest] = N_i_j_dict[key] / N_i
                    N_ij[source, dest] = N_i_j_dict[key]
                    R_i[source] = avg_t[source,0]
                    if dest > n-2:
                        continue
                    K_hat[source, dest] = K[source, dest]
                
                    
        elif milestone.bd:
            bd_results_file = os.path.join(
                me.project.rootdir, milestone.directory, "results.xml")
            
            if not os.path.exists(bd_results_file):
                bd_directory_list_glob = os.path.join(
                    me.project.rootdir, milestone.directory, 
                    "first_hitting_point_distribution", "lig*/")
                bd_directory_list = glob.glob(bd_directory_list_glob)
                combine_fhpd_results(me, milestone, bd_directory_list)
            
            if os.path.exists(bd_results_file):
                transition_probabilities = \
                    browndye_parse_bd_milestone_results(bd_results_file)
                for key in transition_probabilities:
                    value = transition_probabilities[key]
                    if key == "escaped":
                        pass
                    else:
                        K[n-1, n-2] = value
                        bd_milestone_prob = value
            else:
                print("bd results file not found. K-on cannot be computed.")
    
    Q = np.zeros((n,n), dtype=np.double)
    Q_hat = np.zeros((n-1,n-1), dtype=np.double)
    
    for i in range(n-1):
        for j in range(n-1):
            if i == j:
                Q_hat[i,j] = -1.0 / avg_t[i]
                
            else:
                Q_hat[i,j] = K_hat[i,j] / avg_t[i]
            
            Q[i,j] = Q_hat[i,j]
    
    Q[n-2,n-1] = -(Q[n-2,n-3] + Q[n-2,n-2])
    Q[n-1,n-1] = -1.0e-3
    Q[n-1,n-2] = 1.0e-3
    R_i[n-1] = 1e-3
    #K[n-1, n-2] = 1.0
    N_ij[n-1, n-2] = 1e3
                
    return n, K, K_hat, Q_hat, Q, N_ij, R_i, bd_milestone_prob

# function to compute quantities from matrix
def compute_kinetics(n, Q_hat, K, _bd_sample_from_normal=False):
    """
    
    """
    """
    q0 = np.zeros((n-1,1), dtype=np.double)
    q0[0,0] = 1.0
    I = np.matrix(np.identity(n-1))
    aux = np.linalg.solve(I - K_hat, avg_t)
    mfpt = q0.T.dot(aux)[0,0]
    k_off = 1.0e12 / mfpt
    """
    neg_1 = -np.ones((n-1, 1))
    mfpts = np.linalg.solve(Q_hat, neg_1)
    mfpt = mfpts[0,0]
    k_off = 1.0e12 / mfpt
    # k-on
    
    results_filename = os.path.join(
                me.project.rootdir, "b_surface", "results.xml")
    
    k_on = 0
    if os.path.exists(results_filename):
        k_on0, k_on_error, reaction_probability, \
                    reaction_probability_error = \
                    browndye_run_compute_rate_constant(os.path.join(
                        "", "compute_rate_constant"), results_filename, 
                        sample_error_from_normal=_bd_sample_from_normal)
        
        source_vec = np.zeros((n,1))
        source_vec[n-1,0] = k_on0
        """
        K_trap = K[:,:]
        K_trap[0,:] = 0.0
        K_trap[0,0] = 1.0
        K_trap[n-1, :] = 0.0
        K_trap[n-1, n-1] = 1.0
        """
        K_inf = np.linalg.matrix_power(K, MATRIX_EXPONENTIAL)
        end_k_ons = np.dot(K_inf.T, source_vec)
        k_on = end_k_ons[0,0]
        
    else:
        print("No BD files, skipping k-on calculations")
    
    K_stat = K[:,:]
    K_stat[0,:] = 0.0
    K_stat[0,1] = 1.0
    K_stat[n-1, :] = 0.0
    K_stat[n-1, n-2] = 1.0
    
    eigvals, eigvecs = la.eig(K_stat.T)
    closest_eigval = -1e9
    closest_eigvec = None
    for i, eigval, in enumerate(eigvals):
        if (eigval - 1.0)**2 < (closest_eigval - 1.0)**2:
            closest_eigval = eigval
            closest_eigvec = eigvecs[:,i]
    
    p_i = closest_eigvec / np.sum(closest_eigvec)
    
    p_i = p_i @ np.linalg.matrix_power(K_stat, 99) ## !!
    p_i = p_i / np.sum(p_i)
    
    return mfpt, k_off, k_on, p_i

def monte_carlo_milestoning_error(N, R, Q, bd_milestone_prob, num =1000, 
                                  skip=100, stride=1):
        """
        
        """
        m = Q.shape[0] #get size of count matrix
        Q_mats = []
        MFPTs_list = []
        k_off_list = []
        running_avg = []
        running_std = []
        k_on_list = []
        k_on_avg = []
        k_on_std = []
        Qnew = Q[:,:]
        p_i_list = []
        
        for counter in range(num * (stride) + skip):
            #if verbose: print("MCMC stepnum: ", counter)
            Qnew = Q[:,:]
            for i in range(m): # rows
                for j in range(m): # columns
                    if i == j: continue
                    if Qnew[i,j] == 0.0: continue
                    if Qnew[j,j] == 0.0: continue
                    Q_gamma = 0
                    delta = Qnew[i,j]
                    while ((delta) >= (Qnew[i,j])):
                        Q_gamma = gamma.rvs(a=N[i,j], scale = 1/R[i],)
                        delta =  Qnew[i,j] - Q_gamma
        
                    log_p_Q_old = N[i,j] * log(Qnew[i,j])  - Qnew[i,j] * R[i] 
                    log_p_Q_new = N[i,j] * log(Qnew[i,j] - delta) - \
                        (Qnew[i,j] - delta) * R[i]     
                    r2 = random.random()  
                    p_acc =  log_p_Q_new - log_p_Q_old
                        
                    if log(r2) <= p_acc: 
                        #log(r) can be directly compared to 
                        # log-likelihood acceptance, p_acc
                        Qnew[i,i] = (Qnew[i,i]) + delta
                        Qnew[i,j] = Qnew[i,j] - delta
                        
            if counter > skip and counter % stride == 0:
                Knew = np.zeros((m,m))
                for i in range(m):
                    for j in range(m):
                        if i == j: continue
                        if i < n-1:
                            Knew[i,j] = -Qnew[i,j] / Qnew[i,i]
                
                Knew[m-1,m-2] = bd_milestone_prob
                Knew[0,:] = 0.0
                Knew[0,0] = 1.0
                                        
                mfpt, k_off, k_on, p_i = compute_kinetics(n, Qnew[:-1,:-1], Knew, 
                                                     True)
                MFPTs_list.append(mfpt)
                
                k_off_list.append(k_off)
                #if self.model.k_on_info:
                #    for key in k_ons:
                #        k_ons_list[key].append(k_ons[key]) 
                k_on_list.append(k_on) 
                p_i_list.append(p_i)
     
            Q = Qnew
        mfpt_error = np.std(MFPTs_list)
        k_off_error = np.std(k_off_list)
        k_on_error = np.std(k_on_list)
        p_i_error = np.zeros(p_i.shape)
        for i in range(p_i_error.shape[0]):
            p_i_val_list = []
            for j in range(len(p_i_list)):
                p_i_val_list.append(p_i_list[j][i])
            p_i_error[i] = np.std(p_i_val_list)
        #if self.model.k_on_info:
        #    for key in k_ons:
        #        self.k_ons_error[key] = np.std(k_ons_list[key])
        return mfpt_error, k_off_error, k_on_error, p_i_error

def print_results(mfpt, mfpt_err, k_off, k_off_err, k_on, k_on_err, 
                  temperature):
    """
    Print all results of the analysis calculation
    """
    print("k_off (1/s):", pretty_string_value_error(k_off, 
                                                     k_off_err))
    """
    print("k_ons (1/s * 1/M):")
    for key in self.k_ons:
        if key in self.k_ons_error:
            print("  k_on to state", key, ":", pretty_string_value_error(
                float(self.k_ons[key]), float(self.k_ons_error[key])))
        else:
            print("  k_on to state", key, ":", pretty_string_value_error(
                float(self.k_ons[key]), None))
    """
    
    print("Mean first passage time (s):")
    print(pretty_string_value_error(float(mfpt*1.0e-12), mfpt_err*1.0e-12))
    
    print("k_on (s^-1 M^-1):", pretty_string_value_error(k_on, 
                                                     k_on_err))
    diss_constant = k_off / k_on
    delta_G = GAS_CONSTANT*temperature*log(diss_constant)
    diss_constant_err = diss_constant * quadriture(
        k_on_err/k_on, k_off_err/k_off)
    delta_G_err = diss_constant_err*GAS_CONSTANT*temperature/diss_constant
    print("Dissociation constant (M):", 
          pretty_string_value_error(diss_constant, diss_constant_err))
    print("\u0394G (kcal/mol):", 
          pretty_string_value_error(delta_G, delta_G_err))
    return

def save_plots(image_directory, p_i, p_i_error, temperature):
    """
    
    """
    milestone_indices = np.zeros(p_i.shape[0], dtype=np.int8)
    for i in range(p_i.shape[0]):
        milestone_indices[i] = i
    
    # save p_i
    pi_fig, ax = plt.subplots()
    #ax.plot(milestone_indices, p_i, linestyle='-', 
    #        marker="o", markersize = 1)
    plt.errorbar(milestone_indices, p_i, yerr=p_i_error, ecolor="k", capsize=2)
    plt.ylabel("p_i")
    plt.xlabel("milestones")
    pi_fig.savefig(os.path.join(image_directory, "p_i.png"))
    
    free_energy_profile = np.zeros(p_i.shape)
    free_energy_profile_err = np.zeros(p_i.shape)
    highest_p_i = max(p_i)
    for i, p_i_val in enumerate(p_i):
        free_energy = -GAS_CONSTANT*temperature*log(p_i_val / highest_p_i)
        free_energy_err = -GAS_CONSTANT*temperature*highest_p_i*p_i_error[i]\
            /p_i_val
        free_energy_profile[i] = free_energy
        free_energy_profile_err[i] = free_energy_err
    

    # save free energy profile
    pi_fig, ax = plt.subplots()
    #ax.plot(milestone_indices, free_energy_profile, linestyle='-', 
    #        marker="o", markersize = 1)
    plt.errorbar(milestone_indices, free_energy_profile, 
                 yerr=free_energy_profile_err, ecolor="k", capsize=2)
    plt.ylabel("\u0394G(milestone) (kcal/mol)")
    plt.xlabel("milestones")
    pi_fig.savefig(os.path.join(image_directory, "free_energy_profile.png"))
    
    return

if __name__ == "__main__":
    picklename = sys.argv[1]
    assert os.path.exists(picklename), "SEEKR pickle not found in provided "\
            "directory. Are you sure this is a SEEKR calculation directory?"
    
    
    print("Loading SEEKR calculation.")
    me = seekr.openSeekrCalc(picklename)
    image_directory = os.path.join(me.project.rootdir, DEFAULT_IMAGE_DIR)
    if not os.path.exists(image_directory):
        os.mkdir(image_directory)
    n, K, K_hat, Q_hat, Q, N_ij, R_i, bd_milestone_prob = make_matrices(me)
    mfpt, k_off, k_on, p_i = compute_kinetics(n, Q_hat, K)
    mfpt_error, k_off_error, k_on_error, p_i_error = \
        monte_carlo_milestoning_error(N_ij, R_i, Q, bd_milestone_prob, 
                                      num=1000, skip=100, stride=1)
    
    print_results(mfpt, mfpt_error, k_off, k_off_error, k_on, k_on_error, 
                  me.master_temperature)
    print("All plots being saved to:", image_directory)
    save_plots(image_directory, p_i, p_i_error, me.master_temperature)