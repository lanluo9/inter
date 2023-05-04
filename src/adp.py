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
# def load_trace_trial_data(dir_path, vis_filter=True): # deprecated
#     """
#     load trace by trial .mat and stim info data from directory.

#     args:
#         dir_path (raw string): directory path, compatible with windows path containing '\\'
#         vis_filter (boolean): whether to filter by visually driven cells

#     returns:
#         trace_by_trial (ndarray): neural activity of shape: ncell x ntrial x nframe
#         stim_id (ndarray): stim identity for each trial
#     """
#     file_path = dir_path.replace("\\", "/")
#     data = sio.loadmat(file_path + "/trace_trial_stim.mat")
#     # print(data.keys())

#     try:
#         stim_seq = data["stim_seq"]
#         stim_id = [i[0] for i in stim_seq]  # flatten list
#         if stim_seq.shape[1] > stim_seq.shape[0]:  # if stim_seq is a column vector
#             stim_id = stim_seq
#     except:
#         stim_ori = data["stim_ori"]
#         isi_nframe = data["isi_nframe"]
#         adapter_contrast = data["adapter_contrast"]
#         stim_id = dict(
#             stim_ori=stim_ori, isi_nframe=isi_nframe, adapter_contrast=adapter_contrast
#         )
#     trace_by_trial = data["trace_by_trial"]

#     if vis_filter:
#         with open(file_path + "/vis_driven.pickle", "rb") as handle:
#             vis = pickle.load(handle)
#             vis_driven = vis["vis_driven"]
#             vis_driven = [v[0] for v in vis_driven]  # flatten list
#         trace_by_trial = trace_by_trial[
#             vis_driven, :, :
#         ]  # only keep trace of vis-driven cells
#         print("alert! only vis driven cells are kept!")

#     ncell = trace_by_trial.shape[0]
#     try:
#         nstim = len(np.unique(stim_id))
#     except:
#         nstim = len(np.unique(stim_ori))
#     ntrial = trace_by_trial.shape[1]
#     nframe = trace_by_trial.shape[2]
#     print(f"ncell: {ncell}, nstim: {nstim}, ntrial: {ntrial}, nframe: {nframe}")
#     # print(ncell, nstim, ntrial, nframe)

#     return (
#         stim_id,
#         trace_by_trial,
#     )  # ncell, nstim, ntrial, nframe


def load_resp_trial(dir_path, vis_filter=False):
    """
    load resp_cell_trial and stim_info data from directory
    for grat_8ori_3isi paradigm or mix14

    args:
        dir_path (raw string): directory path, compatible with windows path containing '\\'
        vis_filter (boolean): whether to filter by visually driven cells

    returns:
        R1/R2_cell_trial (ndarray): neural activity of shape: ncell x ntrial, avg over resp time window from matlab
        stim_id (ndarray): stim2 identity, isi, stim1 contrast for each trial
    """
    file_path = dir_path.replace("\\", "/")
    data = sio.loadmat(file_path + "/trace_trial_stim.mat")
    # print(data.keys())

    try:
        stim_seq = data["stim_seq"]
        stim_id = [i[0] for i in stim_seq]  # flatten list
        if stim_seq.shape[1] > stim_seq.shape[0]:  # if stim_seq is a column vector
            stim_id = stim_seq
    except:
        stim_ori = data["stim_ori"]
        isi_nframe = data["isi_nframe"]
        adapter_contrast = data["adapter_contrast"]
        stim_id = dict(
            stim_ori=stim_ori, isi_nframe=isi_nframe, adapter_contrast=adapter_contrast
        )
    R1_cell_trial = data["R1_cell_trial"]  # ncell x ntrial, actually not a trace
    R2_cell_trial = data["R2_cell_trial"]

    if vis_filter:
        with open(file_path + "/vis_driven.pickle", "rb") as handle:
            vis = pickle.load(handle)
            vis_driven = vis["vis_driven"]
            vis_driven = [v[0] for v in vis_driven]  # flatten list
        R1_cell_trial = R1_cell_trial[
            vis_driven, :, :
        ]  # only keep trace of vis-driven cells
        R2_cell_trial = R2_cell_trial[
            vis_driven, :, :
        ]  # only keep trace of vis-driven cells
        print("alert! only vis driven cells are kept!")

    ncell = R1_cell_trial.shape[0]
    try:
        nstim = len(np.unique(stim_id["stim_ori"]))  # for grat_8ori_3isi
    except:
        nstim = len(np.unique(stim_id))  # for mix14
    ntrial = R1_cell_trial.shape[1]
    print(f"ncell: {ncell}, nstim: {nstim}, ntrial: {ntrial}")

    return (stim_id, R1_cell_trial, R2_cell_trial)


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
