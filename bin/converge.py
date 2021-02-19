"""
 NOTE TO SELF: finish this if convergence analysis ever necessary...
"""

import os
import argparse
import functools
from collections import defaultdict
from itertools import islice

import numpy as np
import matplotlib.pyplot as plt

import seekr
import analyze

DEFAULT_IMAGE_DIR = "images_and_plots/"
N_IJ_DIR = "N_ij/"
R_I_DIR = "R_i/"
K_CELL_DIR = "k_cell/"
P_EQUIL_DIR = "p_equil/"
MIN_PLOT_NORM = 1e-8

def analyze_kinetics(me, max_step_list):
    
    n, K, K_hat, Q_hat, Q, N_ij, R_i, bd_milestone_prob = \
        analyze.make_matrices(me)
    mfpt, k_off, k_on, p_i = analyze.compute_kinetics(n, Q_hat, K)
    return k_on, k_off, N_ij, R_i

def check_milestone_convergence(me, max_steps):
    conv_stride = max_steps // 100
        
    conv_intervals = np.arange(conv_stride, max_steps, conv_stride)
    max_step_list = np.zeros(len(me.anchors))
    k_off_conv = np.zeros(len(conv_intervals))
    k_on_conv = np.zeros(len(conv_intervals))
    N_ij_conv = defaultdict(functools.partial(np.zeros, len(conv_intervals)))
    R_i_conv = defaultdict(functools.partial(np.zeros, len(conv_intervals)))
    
    for interval_index in range(len(conv_intervals)):
        max_step_list[:] = conv_intervals[interval_index]
        k_on, k_off, N_ij, R_i= analyze_kinetics(me, max_step_list)
        k_on_conv[interval_index] = k_on
        k_off_conv[interval_index] = k_off
        for N_ij_key in N_ij:
            N_ij_conv[N_ij_key][interval_index] = N_ij[N_ij_key] / interval_index
        for R_i_key in R_i:
            R_i_conv[R_i_key][interval_index] = R_i[R_i_key] / interval_index
        
    return k_on_conv, k_off_conv, N_ij_conv, R_i_conv

def plot_rate_const_conv(conv_values, conv_intervals, label="k off (s^-1)"):
    """
    Plot convergence of off/on rate as a function of simulation time

    Parameters
    ----------
    conv_values : list
        list of calculated off rates for each convergence interval
    conv_interval : list
        list of convergence interval step numbers for which samples 
        are taken
    label : str
        The label to give this plot.

    Returns
    -------
    fig : matplotlib figure
        matplotlib figure plotting N convergence for each milestone
    ax : object
        matplotlib Axes object

    """
    fig, ax = plt.subplots()
    ax.plot(np.multiply(conv_intervals,2e-6), conv_values, linestyle='-', 
            marker="o", markersize = 1)
    plt.ylabel(label)
    plt.xlabel("time (ns)")
    return fig, ax

def plot_dict_conv(conv_dict, conv_intervals, label_base="N", skip_null=True):
    """
    Plot convergence of off/on rate as a function of simulation time

    Parameters
    ----------
    conv_dict : dict
        dict of lists of calculated off rates for each convergence 
        interval
    conv_interval : list
        list of convergence interval step numbers for which samples 
        are taken
    label_base : str
        The base of the label to give this plot, though the dictionary
        keys will be appended to the label.
    skip_null : bool
        If true, than empty convergence lists will be omitted from
        any plots.
    
    Returns
    -------
    fig : matplotlib figure
        matplotlib figure plotting N convergence for each milestone
    ax : object
        matplotlib Axes object

    """
    fig_list = []
    ax_list = []
    title_list = []
    for key in conv_dict:
        conv_values = conv_dict[key]
        if skip_null:
            if np.linalg.norm(conv_values) < MIN_PLOT_NORM:
                continue
        if isinstance(key, tuple):
            label = "(" + label_base + " " + "->".join(map(str, key)) + ")"
            title = label_base + "_" + "_".join(map(str, key))
        elif isinstance(key, int):
            label = "(" + label_base + " " + str(key) + ")"
            title = label_base + "_" + str(key)
        else:
            raise Exception("key type not implemented.")
        
        fig, ax = plt.subplots()
        ax.plot(np.multiply(conv_intervals,2e-6), conv_values, linestyle='-', 
                marker="o", markersize = 1)
        plt.ylabel(label)
        plt.xlabel("time (ns)")
        fig_list.append(fig)
        ax_list.append(ax)
        title_list.append(title)
    return fig_list, ax_list, title_list

def _make_windows(seq, n):
    """
    Returns a sliding window (of width n) over data from the iterable
       s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...                   
    """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result
        
def _calc_window_rmsd(conv_values):
    """
    
    """
    average = np.mean(conv_values)
    test = conv_values - average
    new_test = np.delete(test, 0)
    RMSD = np.sqrt(np.sum(new_test)**2/(len(conv_values) -1))
    return RMSD, average

def _find_conv_min(conv_values, conv_intervals, window_size, cutoff, 
                   conv_windows):
    """
    
    """
    conv_times = defaultdict(float)
    for key in conv_values:
        if np.sum(conv_values[key][:]) != 0:
            conv_times[key] = np.nan
            rmsd_list = []
            windows = _make_windows(conv_values[key][:], window_size)
            index = 0
            conv_counter = 0
            for w in windows:
                RMSD, window_average = _calc_window_rmsd(w)
                if RMSD <= (cutoff * window_average):
                    conv_counter += 1
                    if conv_counter == conv_windows:
                        max_int = index + window_size
                        min_time = conv_intervals[max_int]
                        conv_times[key] = min_time
                        break
                else: conv_counter = 0
                index += 1
            
            """    
            if np.isnan(conv_times[key]): print(
                "Entry %s did not meet minimum convergence criteria of "\
                "%i windows" %(str(key), conv_windows))
            """
    return conv_times

def calc_RMSD_conv(model, N_conv, R_conv, conv_intervals, window_size=30, 
                   cutoff=0.01, conv_windows=20):
    """ Estimate the convergence of sampling for each milestone using 
    a sliding window RMSD cutoff. Milestones are considered converged 
    when the values of BOTH N and R remain below the cutoff threshold 
    for a specified number of windows.

    The cutoff is a percentage of the magnutude of the corresponding 
    value (e.g. less than 5% of the magnutude of N).

    Parameters
    -----------
    model : object
        milestoning model object containing all milestone information 
        and transition statistics
    N_conv: list
        list of transition count matrix N for each convergence interval
    R_conv : list
        list of transition time matrix R for each convergence interval 
    conv_intervals : list
        list of stepnumbers for which convergence samples were 
        extracted and calculated
    window size : int
        number of convergence samples used in a window RMSD 
        calculation
    cutoff : float, Default 0.01
        RMSD must be less than the cutoff times the magnitude of the 
        quantity (N or R) to be considered converged
    conv_windows : int, Default 20
        number of consecutive windows the RMSD must remain below the 
        cutoff to be considered converged

    Return
    -------
    min_anchor_times : list
        list of times (in stepnumbers) for each voronoi cell where 
        the convergence criteria are satisfied
    """

    min_anchor_times = np.zeros((len(conv_intervals)))
    #print("Calculating N/T convergence")
    n_conv_times = _find_conv_min(N_conv, conv_intervals, window_size, cutoff,
                                  conv_windows)
    #print("Calculating R/T convergence")
    r_conv_times = _find_conv_min(R_conv, conv_intervals, window_size, cutoff,
                                  conv_windows)
    
    for anchor in model.anchors:
        n_steps = 0
        r_steps = 0
        milestones = []
        for milestone in anchor.milestones:
            milestones.append(milestone.index)
            
        for key in n_conv_times:
            if key[0] in milestones and key[1] in milestones:
                if np.isnan(n_conv_times[key]):
                    n_steps = None
                    break
                if n_conv_times[key] > n_steps:
                    n_steps = n_conv_times[key]
        
        for key in r_conv_times:
            if key in milestones:
                if np.isnan(r_conv_times[key]):
                    r_steps = None
                    break
                if r_conv_times[key] > r_steps:
                    r_steps = r_conv_times[key]
        
        if r_steps is None or n_steps is None:
            print("Anchor %i is not converged." % anchor.index)
        else:
            print("Anchor %i converged after time %f ns" \
                  % (anchor.index, 2e-6 * max(n_steps, r_steps)))
    
    """
    for anchor in model.anchors:
        n_steps = 0
        if np.isnan(r_conv_times[anchor.index]):
            continue
        r_steps = max(r_conv_times[anchor.index])
        for milestone in anchor.milestones:
            if np.isnan(n_conv_times[int(milestone.id)][:]).any():
                continue
            else:
                if max(n_conv_times[int(milestone.id)][:]) > n_steps:               
                    n_steps = max(n_conv_times[int(milestone.id)][:])
        min_anchor_times[int(anchor.index)] = max(n_steps, r_steps) 
        print("anchor", anchor.index, min_anchor_times[int(anchor.index)])
    """

    #return min_anchor_times
    return

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        "input_file", metavar="INPUT_FILE", type=str, 
        help="name of input file for OpenMMVT calculation. This would be the "\
        "XML file generated in the prepare stage.")
    argparser.add_argument(
        "max_steps", metavar="MAX_STEPS", type=int, 
        help="The number of steps to analyze per anchor.")
    
    
    argparser.add_argument(
        "-r", "--RMSD_convergence", dest="RMSD_convergence", 
        default=False, help="Whether to compute RMSD convergence "\
        "in order to determine whether any of the anchors are converged.", 
        action="store_true")
    argparser.add_argument(
        "-w", "--window_size", dest="window_size", default=30, type=int, 
        help="")
    argparser.add_argument(
        "-c", "--cutoff", dest="cutoff", default=0.01, type=float, 
        help="")
    argparser.add_argument(
        "-o", "--conv_windows", dest="conv_windows", default=20, type=int, 
        help="")
    
    args = argparser.parse_args() # parse the args into a dictionary
    args = vars(args)
    xmlfile = args["input_file"]
    max_steps = args["max_steps"]
    rmsd_convergence = args["RMSD_convergence"]
    window_size = args["window_size"]
    cutoff = args["cutoff"]
    conv_windows = args["conv_windows"]
    
    print("Loading SEEKR calculation.")
    me = seekr.openSeekrCalc(xmlfile)
    
    image_directory = os.path.join(me.project.rootdir, DEFAULT_IMAGE_DIR)
    if not os.path.exists(image_directory):
        os.mkdir(image_directory)
    
    if not os.path.exists(os.path.join(image_directory, N_IJ_DIR)):
        os.mkdir(os.path.join(image_directory, N_IJ_DIR))
        
    if not os.path.exists(os.path.join(image_directory, R_I_DIR)):
        os.mkdir(os.path.join(image_directory, R_I_DIR))
        
    if not os.path.exists(os.path.join(image_directory, K_CELL_DIR)):
        os.mkdir(os.path.join(image_directory, K_CELL_DIR))
        
    if not os.path.exists(os.path.join(image_directory, P_EQUIL_DIR)):
        os.mkdir(os.path.join(image_directory, P_EQUIL_DIR))
        
    k_on_conv, k_off_conv, N_ij_conv, R_i_conv, k_alpha_beta_conv, \
        pi_alpha_conv, conv_intervals = check_milestone_convergence(
            me, max_steps)
    
    k_off_fig, ax = plot_rate_const_conv(k_off_conv, conv_intervals)
    k_on_fig, ax = plot_rate_const_conv(k_on_conv, conv_intervals, 
                                        label="k on (s^-1 M^-1)")
    N_ij_fig_list, ax, N_ij_title_list = plot_dict_conv(
        N_ij_conv, conv_intervals, label_base="N")
    R_i_fig_list, ax, R_i_title_list = plot_dict_conv(
        R_i_conv, conv_intervals, label_base="R")
    k_alpha_beta_fig_list, ax, k_alpha_beta_title_list = plot_dict_conv(
        k_alpha_beta_conv, conv_intervals, label_base="k")
    pi_alpha_fig_list, ax, pi_alpha_title_list = plot_dict_conv(
        pi_alpha_conv, conv_intervals, label_base="pi")
    
    k_off_fig.savefig(os.path.join(image_directory, "k_off_convergence.png"))
    k_on_fig.savefig(os.path.join(image_directory, "k_on_convergence.png"))
    for i, title in enumerate(N_ij_title_list):
        N_ij_fig = N_ij_fig_list[i]
        N_ij_fig.savefig(os.path.join(image_directory, N_IJ_DIR, 
                                      title+".png"))
    for i, title in enumerate(R_i_title_list):
        R_i_fig = R_i_fig_list[i]
        R_i_fig.savefig(os.path.join(image_directory, R_I_DIR, 
                                      title+".png"))
    for i, title in enumerate(k_alpha_beta_title_list):
        k_alpha_beta_fig = k_alpha_beta_fig_list[i]
        k_alpha_beta_fig.savefig(os.path.join(image_directory, K_CELL_DIR, 
                                      title+".png"))
    for i, title in enumerate(pi_alpha_title_list):
        pi_alpha_fig = pi_alpha_fig_list[i]
        pi_alpha_fig.savefig(os.path.join(image_directory, P_EQUIL_DIR, 
                                      title+".png"))
    
    print("All plots have been saved to:", image_directory)
    
    if rmsd_convergence:
        print("RMSD convergence report:")
        min_anchor_times = calc_RMSD_conv(
            me, N_ij_conv, R_i_conv, conv_intervals, window_size, cutoff, 
            conv_windows)