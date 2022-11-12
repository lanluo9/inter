# module docs
__doc__ = """
This module contains functions for the analysis of adaptation.
"""

# import

import numpy as np
import pandas as pd
import seaborn as sns

# sns.set_theme(style="whitegrid")
sns.set_context("talk", font_scale=0.8)

import scipy.io as sio
from scipy import stats
import os
import pickle
from tqdm import tqdm

import matplotlib

matplotlib.rc("xtick", labelsize=16)
matplotlib.rc("ytick", labelsize=16)
import matplotlib.pyplot as plt

plt.rcParams.update({"font.size": 10})

# functions
def load_trace_trial_data(dir_path, vis_filter=True):
    """
    load trace by trial .mat and stim info data from directory.

    args:
        dir_path (raw string): directory path, compatible with windows path containing '\\'
        vis_filter (boolean): whether to filter by visually driven cells

    returns:
        trace_by_trial (ndarray): neural activity of shape: ncell x ntrial x nframe
        stim_id (ndarray): stim identity for each trial
    """
    file_path = dir_path.replace("\\", "/")
    data = sio.loadmat(file_path + "/trace_trial_stim.mat")

    stim_seq = data["stim_seq"]
    stim_id = [i[0] for i in stim_seq]  # flatten list
    trace_by_trial = data["trace_by_trial"]

    if vis_filter:
        with open(file_path + "/vis_driven.pickle", "rb") as handle:
            vis = pickle.load(handle)
            vis_driven = vis["vis_driven"]
            vis_driven = [v[0] for v in vis_driven]  # flatten list
        trace_by_trial = trace_by_trial[
            vis_driven, :, :
        ]  # only keep trace of vis-driven cells

    ncell = trace_by_trial.shape[0]
    nstim = len(np.unique(stim_id))
    ntrial = trace_by_trial.shape[1]
    nframe = trace_by_trial.shape[2]
    print(f"ncell: {ncell}, nstim: {nstim}, ntrial: {ntrial}, nframe: {nframe}")
    # print(ncell, nstim, ntrial, nframe)

    return (
        stim_id,
        trace_by_trial,
    )  # ncell, nstim, ntrial, nframe


def calc_trace_stim(trace_by_trial, stim_id):
    """
    calculate mean trace and mean trace by stim.

    args:
        trace_by_trial (ndarray): neural activity of shape: ncell x ntrial x nframe
        stim_id (ndarray): stim identity for each trial

    returns:
        trace_mean (ndarray): mean neural activity of shape: ncell x nframe
        trace_by_stim (ndarray): mean neural activity of shape: nstim x ncell x nframe
    """

    # trace over all stims
    trace_cell_avg = np.mean(np.mean(trace_by_trial, axis=0), axis=0)
    trace_cell_sem = np.std(np.mean(trace_by_trial, axis=0), axis=0) / np.sqrt(
        trace_by_trial.shape[0]
    )

    # trace for each stim
    trace_stim_avg = []
    for i in np.unique(stim_id):
        trace_istim_avg = np.mean(trace_by_trial[:, np.where(stim_id == i)[0]], axis=1)
        trace_istim_avg = np.mean(trace_istim_avg, axis=0)
        trace_stim_avg.append(trace_istim_avg)
    print(
        f"trace_cell_avg: {trace_cell_avg.shape}. \
        trace_stim_avg list len: {len(trace_stim_avg)}. \
        trace_stim_avg[0].shape: {trace_stim_avg[0].shape}"
    )

    return trace_cell_avg, trace_cell_sem, trace_stim_avg


def read_csv_by_stim_type():
    """
    read csv file containing metadata of recordings, including stim type 

    returns:
        df (pandas dataframe): dataframe of recordings for each stim type
    """
    # read metadata of segmented sessions
    dir_inter = r"C:\Users\ll357\Documents\inter\data".replace("\\", "/")
    df = pd.read_csv(dir_inter + "/batch_cellpose.csv")
    # only keep segmented data
    df = df[(df.manual_seg == 1) | (df.cellpose_seg == 1)].reset_index(drop=True)
    # separate by stim type
    df_bun = df[df["stim_type"] == "bunny"]
    df_grat = df[(df["stim_type"] == "grating")]
    return df, df_bun, df_grat
