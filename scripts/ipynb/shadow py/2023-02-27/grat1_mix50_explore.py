#!/usr/bin/env python
# coding: utf-8

# # import

# In[506]:


import numpy as np
import pandas as pd
# import scipy.io
# from scipy import stats

import seaborn as sns
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

# from tqdm import tqdm
import os
import pickle

get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')


# In[507]:


local_flag = False
if local_flag:
    repo_dir = r'D:\repo\inter_data\inter'.replace("\\", "/") # under env dimred
else:
    repo_dir = r'C:\Users\ll357\Documents\inter'.replace("\\", "/")
os.chdir(repo_dir)
from src import adp


# # load mix50
# mix50 stim: Z:\All_Staff\home\lan\Mwork\mix50 - bunnytop high lum contrast mix grating and noise\Image

# In[508]:


if local_flag:
    dir_data = r'D:\repo\inter_data\mix50\V1_i1373_220909_cellpose'
else:
    dir_data = r'\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images\mat_inter\mix50_data/'

# Create a dictionary with keys being recording folder name, and value -- array with first element being stimulus identity array, and second element -- trace by trial matrix; 
recordings = {f : [adp.load_trace_trial_data(os.path.join(dir_data, f), vis_filter=True)] for f in os.listdir(dir_data) if os.path.isdir(os.path.join(dir_data, f))}

# Calculate the minimum number of trials for each trial id across recordings and store in a dictionary with trial id as the key, and smallest number of trials as a value
min_trial_num = {}
for recording in recordings.keys():
    si, tbt = recordings[recording][0]
    names, counts = np.unique(si, return_counts=True)
    for i, n in enumerate(names):
        if n not in min_trial_num.keys() or counts[i] < min_trial_num[n]:
            min_trial_num[n] = counts[i]


# In[509]:


#out_tbt is the final merged array
out_tbt = np.array([])
trialids = min_trial_num.keys()

for recording in recordings.keys():
    #for each recording, re-create trace-by-trial matrix with trial ids sorted the same way across recordings
    out_tbt_by_recording = np.array([])
    for trial_id in trialids:
        min_num_trials = min_trial_num[trial_id]
        # print(min_num_trials)
        si, tbt = recordings[recording][0]
        #select the trials where the current trial_id is used and use that to index
        curr_trialid_locs = np.where(si == trial_id)[0][:min_num_trials]
        # print(curr_trialid_locs.shape)
        # print('\n')

        tbt_slice = tbt[:, curr_trialid_locs, :]
        if out_tbt_by_recording.shape[0] == 0:
            out_tbt_by_recording = tbt_slice
        else:
            #stack trace-by-trial matrices along 1st dimension (n trial)
            out_tbt_by_recording = np.hstack((out_tbt_by_recording, tbt_slice))
    #stack trace-by-trial matrices along 0th dimesion (neuron)
    if out_tbt.shape[0] == 0:
        out_tbt = out_tbt_by_recording
    else:
        out_tbt = np.vstack((out_tbt, out_tbt_by_recording))

print(out_tbt.shape)


# In[510]:


stim_id_merged = np.array([])
for trial_id in trialids:
    min_num_trials = min_trial_num[trial_id]
    for i in range(min_num_trials):
        stim_id_merged = np.append(stim_id_merged, int(trial_id))
stim_id_merged = stim_id_merged.astype(int)

# count stim_id_merged elements
np.unique(stim_id_merged, return_counts=True)
stim_id_merged.shape


# In[511]:


trace_by_trial = out_tbt
stim_id = stim_id_merged
trace_avg_cell, trace_cell_sem, trace_stim_avg = adp.calc_trace_stim(trace_by_trial, stim_id)

fig, ax = plt.subplots(figsize=(10, 5))
plt.plot(trace_avg_cell[:50], alpha=0.7, linewidth=3)
ax.axvspan(9,11, alpha=0.2, color='red')
ax.axvspan(13,14, alpha=0.2, color='gray')
ax.axvspan(20,22, alpha=0.2, color='red')
plt.xlim(0, 50);


# In[512]:


stim_type_dict = {'natural': [np.arange(1, 30+1)], 
                  'grat': [np.arange(31, 40+1)], 
                  'noise': [np.arange(41, 50+1)]}
for key in stim_type_dict.keys():
    print(key, np.min(stim_type_dict[key]), '-', np.max(stim_type_dict[key]))


# In[513]:


# get resp_ad and resp_tg [cell x trial]
base1 = trace_by_trial[:,:,0:2+1].mean(axis=2) # avg over time window frames
resp_ad = trace_by_trial[:,:,9:11+1].mean(axis=2)
resp_ad = resp_ad - base1

base2 = trace_by_trial[:,:,13:14+1].mean(axis=2)
resp_tg = trace_by_trial[:,:,20:22+1].mean(axis=2)
resp_tg = resp_tg - base2

resp_ad.shape, resp_tg.shape


# # vis, img driven, pref filter
# merged mix50 data only support 'vis' filter for now

# In[514]:


img_driven_merge = np.array([])

for subdir in os.listdir(dir_data):
    data_dir = os.path.join(dir_data, subdir)
    print(data_dir)

    with open(data_dir + "/vis_driven.pickle", "rb") as handle:
        vis = pickle.load(handle)
    vis_driven = vis["vis_driven"]
    vis_driven = [v[0] for v in vis_driven]

    img_driven = vis['img_driven'][vis_driven, :]
    print(img_driven.shape)
    img_driven_merge = img_driven if img_driven_merge.size == 0 else np.concatenate((img_driven_merge, img_driven))

img_driven = img_driven_merge

img_driven = np.ones_like(img_driven)
print('turn off img_driven for now')
img_driven.shape


# In[515]:


cell_filter_name = 'img_driven'
print(cell_filter_name)
if cell_filter_name == 'img_driven': # only use img driven cells
    cell_filter = img_driven
# elif cell_filter_name == 'pref_filter': # only use cells preferring each img
#     pref_filter = np.load(data_dir + "/pref_filter.npy") # load pref_filter.npy, constructed under img driven filter
#     cell_filter = pref_filter.copy().T
elif cell_filter_name == 'vis': # use all vis driven cells
    cell_filter = np.ones_like(img_driven)
else:
    raise ValueError('cell_filter not recognized')

plt.matshow(cell_filter.T)
plt.xticks([]) # no axis ticks
plt.yticks([])
plt.xlabel('cells')
plt.ylabel('stims');

print(cell_filter.shape) # ncell x nstim
print(f'{np.sum(cell_filter.sum(axis=0) == 0)} out of {cell_filter.shape[1]} stims have no cells passing filter')
print(f'{sum(cell_filter.sum(axis=1) == 0)} out of {cell_filter.shape[0]} cells have no stim passing filter')


# In[516]:


def resp_by_stim(resp_ad, resp_tg, cell_filter, stim_id):
    '''
    reshape resp_ad and resp_tg by stim: nstim x [ncell x trial_rep]
    enable stim_id to be modified, so we can bin spatial freq of gratings
    '''

    resp_ad_by_stim = []
    resp_tg_by_stim = []

    for istim in np.unique(np.sort(stim_id)):
        istim_id = np.where(stim_id == istim)[0] # which trial has this stim
        
        resp_ad_istim = resp_ad[:, istim_id]
        resp_ad_istim[cell_filter[:, istim-1] == False, :] = np.nan # stim id starts from 1, but index starts from 0
        # TODO: comment line above if using mix50 merged data
        resp_ad_by_stim.append(resp_ad_istim)
        
        resp_tg_istim = resp_tg[:, istim_id]
        resp_tg_istim[cell_filter[:, istim-1] == False, :] = np.nan
        resp_tg_by_stim.append(resp_tg_istim)

    return resp_ad_by_stim, resp_tg_by_stim

resp_ad_by_stim, resp_tg_by_stim = resp_by_stim(resp_ad, resp_tg, cell_filter, stim_id)
# len(resp_ad_by_stim), [i.shape for i in resp_ad_by_stim[:5]]
# resp_ad_by_stim[0].shape, len(resp_ad_by_stim) # ncell x ntrial_rep x nstim
resp_ad_stim_cell = np.array([np.nanmean(i, axis=1) for i in resp_ad_by_stim]) # [ncell x nstim], mean over trial reps
resp_ad_stim_cell.shape # nstim x ncell

resp_ad_cell_stim = resp_ad_stim_cell.T
resp_ad_cell_stim.shape # ncell x nstim, aka nsample x nfeature
# resp_ad_cell_stim[np.isnan(resp_ad_cell_stim)] = 0 # substitute nan with 0, to accomodate hierarchical clustering

plt.matshow(resp_ad_cell_stim.T)
plt.xlabel('cells')
plt.ylabel('stimuli');


# # yuansi regression adaptation estimator

# ## assumption 1: gaussian
# R1 and R2 are approx gaussian distributed, given a certain stim

# In[549]:


def melt_data(resp_mat):
    '''
    melt data into long format
    resp_mat # (ncell x ntrial_rep) x nstim
    '''
    cell_id_concat = []
    stim_id_concat = []
    resp_concat = []

    for istim in range(len(resp_mat)):
        resp_stim = resp_mat[istim] # ncell x ntrial_rep
        for icell in range(resp_stim.shape[0]):
            resp_cell = resp_stim[icell, :] # ntrial_rep
            cell_id_rep = np.repeat(icell, len(resp_cell)) # ntrial_rep
            stim_id_rep = np.repeat(istim, len(resp_cell)) # ntrial_rep
            cell_id_concat.append(cell_id_rep)
            stim_id_concat.append(stim_id_rep)
            resp_concat.append(resp_cell)

    return cell_id_concat, stim_id_concat, resp_concat

cell_id_concat, stim_id_concat, resp_concat = melt_data(resp_ad_by_stim)
df_R1 = pd.DataFrame({'cell_id': np.concatenate(cell_id_concat), 'stim_id': np.concatenate(stim_id_concat), 'R1': np.concatenate(resp_concat)})

cell_id_concat, stim_id_concat, resp_concat = melt_data(resp_tg_by_stim)
df_R2 = pd.DataFrame({'cell_id': np.concatenate(cell_id_concat), 'stim_id': np.concatenate(stim_id_concat), 'R2': np.concatenate(resp_concat)})

assert sum(df_R1.cell_id == df_R2.cell_id) / len(df_R1.cell_id) == 1
assert sum(df_R1.stim_id == df_R2.stim_id) / len(df_R1.stim_id) == 1

df_resp = df_R1.copy()
df_resp['R2'] = df_R2['R2']
df_resp


# In[655]:


# on population level, given a stim, is R1 and R2 gaussian distributed? --> no, also not lognormal after right shift, also not normal after sqrt
# on cell level, given a stim, is R1 and R2 gaussian distributed? --> R2 might be gauss for some cells (like cell 42)? R1 not

df_resp_now = df_resp.copy()
stim_id_now = 0 # TODO: loop over stim_id
df_resp_now = df_resp_now[df_resp_now.stim_id == stim_id_now]
cell_id_now = 42 # TODO: loop over cell_id
df_resp_now = df_resp_now[df_resp_now.cell_id == cell_id_now]
R1_stim = df_resp_now['R1'].values
R2_stim = df_resp_now['R2'].values

# # log resp to check if lognormal
# R1_stim = R1_stim - np.amin(R1_stim) + 1e-3 # shift to positive to log
# R2_stim = R2_stim - np.amin(R2_stim) + 1e-3
# R1_stim = np.log(R1_stim)
# R2_stim = np.log(R2_stim)

# # sqrt resp to check if normal
# R1_stim = np.sqrt(R1_stim)
# R2_stim = np.sqrt(R2_stim)

plt.hist(R1_stim, bins=50, alpha=0.5)
plt.hist(R2_stim, bins=50, alpha=0.5)
plt.legend(['R1', 'R2']);


# In[656]:


import statsmodels.api as sm

fig, ax = plt.subplots(1, 2, figsize=(10, 5))
sm.qqplot(R1_stim, line='45', ax=ax[0]) # Q-Q plot with 45-degree line 
sm.qqplot(R2_stim, line='45', ax=ax[1])
plt.show()


# In[657]:


# https://www.statology.org/normality-test-python/
# https://towardsdatascience.com/normality-tests-in-python-31e04aa4f411

from scipy.stats import normaltest
print(normaltest(R1_stim, axis=0, nan_policy='propagate'))
print(normaltest(R2_stim, axis=0, nan_policy='propagate')) # Test whether a sample differs from a normal distribution

from scipy.stats import shapiro 
print(shapiro(R1_stim))
print(shapiro(R2_stim)) # Shapiro-Wilk test

from scipy.stats import kstest # c:\Users\ll357\Anaconda3\lib\site-packages\scipy\stats\morestats.py:1760: UserWarning: p-value may not be accurate for N > 5000.
print(kstest(R1_stim, 'norm')) # Kolmogorov-Smirnov test
print(kstest(R2_stim, 'norm'))


# ## assumption 2: noise independence
# R1 and R2 are independent in terms of their variance: there is no correlation between variance of R1 and R2

# In[659]:


# on population level, given a stim, is R1 and R2 var independent? --> scatter looks correlated: even if indep, will have multicollinearity

R1_stim_var_list = []
R2_stim_var_list = []

for istim in (df_resp.stim_id.unique()):
    stim_id_now = istim
    df_resp_now = df_resp[df_resp.stim_id == stim_id_now]
    cell_id_now = 42
    df_resp_now = df_resp_now[df_resp_now.cell_id == cell_id_now]

    R1_stim = df_resp_now['R1'].values
    R2_stim = df_resp_now['R2'].values

    R1_stim_var = np.var(R1_stim)
    R2_stim_var = np.var(R2_stim)
    R1_stim_var_list.append(R1_stim_var)
    R2_stim_var_list.append(R2_stim_var)

R1_stim_var_arr = np.array(R1_stim_var_list)
R2_stim_var_arr = np.array(R2_stim_var_list)
plt.scatter(R1_stim_var_arr, R2_stim_var_arr);


# In[651]:


# on cell level, given a stim, is R1 and R2 var independent? --> 

# ncell = len(df_resp.cell_id.unique())
ncell = 15
# set up fig with subplots, one for each cell
fig, ax = plt.subplots(ncell//5+1, 5, figsize=(20, 5))

for icell in range(ncell):
    df_resp_now = df_resp.copy()
    cell_id_now = icell
    df_resp_now = df_resp_now[df_resp_now.cell_id == cell_id_now]

    R1_stim_var_list = []
    R2_stim_var_list = []
    for istim in (df_resp_now.stim_id.unique()):
        stim_id_now = istim
        df_resp_now = df_resp_now[df_resp_now.stim_id == stim_id_now]

        R1_stim = df_resp_now['R1'].values
        R2_stim = df_resp_now['R2'].values
        print(R1_stim.shape, R2_stim.shape)

        R1_stim_var = np.var(R1_stim)
        R2_stim_var = np.var(R2_stim)
        R1_stim_var_list.append(R1_stim_var)
        R2_stim_var_list.append(R2_stim_var)

    R1_stim_var_arr = np.array(R1_stim_var_list)
    R2_stim_var_arr = np.array(R2_stim_var_list)
    ax[icell//5, icell%5].scatter(R1_stim_var_arr, R2_stim_var_arr)


# ## assumption 3: linear
# linear relationship: stim feature -> R1_mean   
# linear relationship: stim feature -> (R2_mean / R1_mean)

# In[625]:


# TODO: quantify stim feature, human defined or ML defined - see workflowy
# TODO: then scatter w R1_mean or R2_mean/R1_mean


# # hierarchical clustering
# find neurons w similar resp profile or adp profile  
# compare neuron clustering based on resp vs adp

# In[492]:


# source: https://scikit-learn.org/stable/auto_examples/cluster/plot_agglomerative_dendrogram.html
from scipy.cluster.hierarchy import dendrogram
from sklearn.cluster import AgglomerativeClustering


def plot_dendrogram(model, **kwargs):
    # Create linkage matrix and then plot the dendrogram

    # create the counts of samples under each node
    counts = np.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)
    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [model.children_, model.distances_, counts]
    ).astype(float)

    # Plot the corresponding dendrogram
    dendrogram(linkage_matrix, **kwargs)

    return linkage_matrix


# ## clustering by response

# ### grating tuning cluster

# In[493]:


def cluster_stim_resp(distance_threshold, resp_stim_type, nlevel_visible, stim_type_str):
    model = AgglomerativeClustering(distance_threshold=distance_threshold, n_clusters=None, compute_distances=False) # set distance_threshold to non zero to get model.labels_ usable - we dont want every cell to be its own cluster
    model = model.fit(resp_stim_type)

    _, axes = plt.subplots(1, 1, figsize=(20, 10))
    _ = plot_dendrogram(model, ax=axes, truncate_mode="level", p=nlevel_visible, count_sort='descending', leaf_rotation=90, leaf_font_size=12) # plot the first p levels of the dendrogram
    plt.axhline(y=distance_threshold, c='gray', alpha=0.5, linestyle='--')

    plt.title(f"Hierarchical Clustering Dendrogram, {stim_type_str}")
    plt.xlabel("Number of points in node (or index of point if no parenthesis).")
    plt.ylabel('Distance');
    # plt.savefig(f'dendrogram_{stim_type_str}_tuning.pdf')

    print(type(model))
    return model


# In[494]:


def norm_gauss_across_stim(resp_mat):
    # resp_mat shape: ncell x nstim
    resp_mat = resp_mat - np.nanmean(resp_mat, axis=1)[:, None] # normalize resp_grat across stims, but this ends up with weird distribution
    resp_mat = resp_mat / np.nanstd(resp_mat, axis=1)[:, None]
    resp_mat[np.isnan(resp_mat)] = 0 # substitute nan with 0, to accomodate hierarchical clustering
    return resp_mat


def norm_maxmin_across_stim(adp_mat):
    # max min normalize resp or adp_mat across stims, mat shape: ncell x nstim
    adp_mat = (adp_mat - np.nanmin(adp_mat, axis=1)[:, None]) / (np.nanmax(adp_mat, axis=1)[:, None] - np.nanmin(adp_mat, axis=1)[:, None])
    adp_mat[np.isnan(adp_mat)] = 0 # substitute nan with 0, to accomodate hierarchical clustering
    return adp_mat


# In[495]:


resp_grat = resp_ad_cell_stim[:, 30:40]
resp_grat = norm_gauss_across_stim(resp_grat)
plt.hist(resp_grat.flatten(), bins=100);
# plt.matshow(resp_grat.T < 0)


# In[498]:


resp_grat = resp_ad_cell_stim[:, 30:40]
resp_grat = norm_maxmin_across_stim(resp_grat)
# plt.matshow(resp_grat.T)
model_grat = cluster_stim_resp(3, resp_grat, 5, 'gratings')


# In[499]:


# # outlier_cells = np.array([356, 231, 164, 347]) # no normalization
# outlier_cells = np.array([175, 157, 94]) # with normalization

# for icell in outlier_cells:
#     plt.plot(resp_grat[icell, :], alpha=0.5, label=str(icell))
# plt.legend()
# plt.xlabel('grating SF')
# plt.ylabel('cell tuning curve');


# In[500]:


# # cluster_size_sort = np.argsort(np.argmax(resp_grat, axis=1))[::-1]
# # cluster_size_sort,
# # resp_grat.shape
# # # cluster_size_sort = np.unique(model_grat.labels_)[cluster_size_sort]

# # sort clusters by first peak in tuning curve
# resp_grat_agg = np.zeros((len(np.unique(model_grat.labels_)), resp_grat.shape[1]))
# for ilabel in np.unique(model_grat.labels_):
#     resp_grat_agg[ilabel, :] = np.nanmean(resp_grat[model_grat.labels_ == ilabel, :], axis=0)
# cluster_peak_sort = np.argmax(resp_grat_agg, axis=1) # for cluster, find max resp stim
# cluster_peak_sort = np.argsort(cluster_peak_sort)
# cluster_peak_sort = np.unique(model_grat.labels_)[cluster_peak_sort]
# # print(cluster_peak_sort)


# In[501]:


def plot_cluster_tuning(model, resp_stim_type, stim_type_str, cluster_sortby='size', sharexy=True, ncol=5, figsize=(15, 15)):

    nsubplot = len(np.unique(model.labels_))
    nrow = nsubplot // ncol + 1 # set 5 columns, n rows

    if sharexy:
        _, ax = plt.subplots(nrow, ncol, figsize=figsize, sharex=True, sharey=True)
    else:
        _, ax = plt.subplots(nrow, ncol, figsize=figsize)

    if cluster_sortby == 'size':
        cluster_size = np.array([len(np.where(model.labels_ == i)[0]) for i in np.unique(model.labels_)])
        cluster_size_sort = np.argsort(cluster_size)[::-1]
        cluster_size_sort = np.unique(model.labels_)[cluster_size_sort] # cluster id sorted by size
        cluster_sort = cluster_size_sort

    elif cluster_sortby == 'tuning peak':
        # sort clusters by first peak in tuning curve
        resp_agg = np.zeros((len(np.unique(model.labels_)), resp_stim_type.shape[1]))
        for ilabel in np.unique(model.labels_):
            resp_agg[ilabel, :] = np.nanmean(resp_stim_type[model.labels_ == ilabel, :], axis=0)
        cluster_peak_sort = np.argmax(resp_agg, axis=1) # for cluster, find max resp stim
        cluster_peak_sort = np.argsort(cluster_peak_sort)
        cluster_peak_sort = np.unique(model.labels_)[cluster_peak_sort]
        cluster_sort = cluster_peak_sort

    for isubplot, ilabel in enumerate(cluster_sort): # for each cluster, subplot sorted by cluster peak tuning
        resp_agg = np.mean(resp_stim_type[model.labels_ == ilabel, :], axis=0) # set nan to 0 before, so should use mean but not median

        if stim_type_str == 'all': # plot diff stim type from same starting point 0
            grat_id = np.arange(30, 40)
            noise_id = np.arange(40, 50)
            nat_id = np.arange(0, 30)
            ax[isubplot // ncol, isubplot % ncol].plot(resp_agg[grat_id], linewidth=4, color='blue', alpha=0.6)
            ax[isubplot // ncol, isubplot % ncol].plot(resp_agg[noise_id], linewidth=4, color='green', alpha=0.6)
            ax[isubplot // ncol, isubplot % ncol].plot(resp_agg[nat_id], linewidth=4, color='orange', alpha=0.6)

            for icell in np.where(model.labels_ == ilabel)[0]:
                ax[isubplot // ncol, isubplot % ncol].plot(resp_stim_type[icell, grat_id], alpha=0.4, linestyle='--', linewidth=1, color='b')
                ax[isubplot // ncol, isubplot % ncol].plot(resp_stim_type[icell, noise_id], alpha=0.4, linestyle='--', linewidth=1, color='g')
                ax[isubplot // ncol, isubplot % ncol].plot(resp_stim_type[icell, nat_id], alpha=0.4, linestyle='--', linewidth=1, color='y')
        
        elif stim_type_str == 'SF':
            grat_id = np.arange(0, 10)
            noise_id = np.arange(10, 20)
            ax[isubplot // ncol, isubplot % ncol].plot(resp_agg[grat_id], linewidth=4, color='blue', alpha=0.6)
            ax[isubplot // ncol, isubplot % ncol].plot(resp_agg[noise_id], linewidth=4, color='green', alpha=0.6)

            for icell in np.where(model.labels_ == ilabel)[0]:
                ax[isubplot // ncol, isubplot % ncol].plot(resp_stim_type[icell, grat_id], alpha=0.4, linestyle='--', linewidth=1, color='b')
                ax[isubplot // ncol, isubplot % ncol].plot(resp_stim_type[icell, noise_id], alpha=0.4, linestyle='--', linewidth=1, color='g')
        else:
            # plot the mean tuning curve across cells in a subplot
            ax[isubplot // ncol, isubplot % ncol].plot(resp_agg, linewidth=4, color='cyan', alpha=0.8)

            # plot the tuning curve of each cell in the cluster
            for icell in np.where(model.labels_ == ilabel)[0]:
                ax[isubplot // ncol, isubplot % ncol].plot(resp_stim_type[icell, :], alpha=0.5, linestyle='--', linewidth=1)

        # add scale bar on the left, if not sharex and sharey
        if not sharexy:
            scale_bar_len = (np.amax(resp_stim_type) - np.amin(resp_stim_type)) * 0.1
            ax[isubplot // ncol, isubplot % ncol].plot([-1, -1], [0, scale_bar_len], linewidth=5, color='gray')

        # add text annotation of cluster id & cell number
        ax[isubplot // ncol, isubplot % ncol].text(0, 1, f'cluster {isubplot}, ncell={np.sum(model.labels_ == ilabel)}', transform=ax[isubplot // ncol, isubplot % ncol].transAxes, fontsize=12, verticalalignment='bottom')

        ax[isubplot // ncol, isubplot % ncol].set_xticks([])
        ax[isubplot // ncol, isubplot % ncol].set_yticks([])
        ax[isubplot // ncol, isubplot % ncol].axis('off')
        # ax[isubplot // ncol, isubplot % ncol].axhline(y=0, color='gray', linewidth=0.5, alpha=0.5) # due to norm, curve min is not 0

    # hide subplots that are not used
    for i in range(nsubplot, nrow * ncol):
        ax[i // ncol, i % ncol].axis('off')

    plt.suptitle(f'Clustered {stim_type_str} tuning curves', fontsize=14, y=0.99)
    plt.tight_layout()
    plt.savefig(f'cluster_{stim_type_str}_tuning.pdf')


# In[502]:


plot_cluster_tuning(model_grat, resp_grat, 'gratings', cluster_sortby='tuning peak', sharexy=False)
# TODO: add error bar for tuning curve to visliz if enough trial rep
# TODO: k means clustering if don't need tree structure or distance 
# TODO: compare clustering by confusion matrix - pretend 1 clustering result is ground truth, 
# and compare the other clustering result to it

# TODO: clustering stability metric - silhouette score, calinski harabasz score, davies bouldin score
# TODO: also used to choose best number of clusters


# ### noise tuning cluster

# In[504]:


resp_noise = resp_ad_cell_stim[:, 40:]
resp_noise = norm_maxmin_across_stim(resp_noise)

model_noise = cluster_stim_resp(3, resp_noise, 5, 'noise')


# In[505]:


plot_cluster_tuning(model_noise, resp_noise, 'noise', cluster_sortby='tuning peak', sharexy=False)


# ### natural image tuning cluster

# In[469]:


resp_nat = resp_ad_cell_stim[:, 0:10]
print('resp nat subsample to 10 stims to match gratings')
resp_nat = norm_maxmin_across_stim(resp_nat)

model_nat = cluster_stim_resp(2, resp_nat, 5, 'natual images')


# In[470]:


plot_cluster_tuning(model_nat, resp_nat, 'natural images', cluster_sortby='tuning peak', sharexy=False)


# ### special cells

# In[472]:


# outlier_cells = np.array([164, ]) # 347

# plt.figure(figsize=(10, 5))
# plt.plot(resp_grat[icell, :], alpha=0.5, label='gratings')
# plt.plot(resp_noise[icell, :], alpha=0.5, label='noise')
# plt.plot(resp_nat[icell, :], alpha=0.5, label='natural images')
    
# plt.legend(frameon=False, bbox_to_anchor=(0.8, 1.05))
# plt.xlabel('stim id')
# plt.ylabel('cell tuning curve');
# plt.gca().spines['top'].set_visible(False)
# plt.gca().spines['right'].set_visible(False)


# In[473]:


# all_cells_exclude_outlier = np.setdiff1d(np.arange(resp_ad_cell_stim.shape[0]), outlier_cells)
# resp_grat_exclude_outlier = resp_grat[all_cells_exclude_outlier, :].mean(axis=0)
# resp_noise_exclude_outlier = resp_noise[all_cells_exclude_outlier, :].mean(axis=0)
# resp_nat_exclude_outlier = resp_nat[all_cells_exclude_outlier, :].mean(axis=0)

# plt.figure(figsize=(10, 5))
# plt.plot(resp_grat_exclude_outlier, alpha=0.5, label='gratings')
# plt.plot(resp_noise_exclude_outlier, alpha=0.5, label='noise')
# plt.plot(resp_nat_exclude_outlier, alpha=0.5, label='natural images')

# plt.legend(frameon=False, bbox_to_anchor=(0.8, 1.05))
# plt.xlabel('stim id')
# plt.ylabel('cell tuning curve');
# plt.gca().spines['top'].set_visible(False)
# plt.gca().spines['right'].set_visible(False)


# ### SF stim tuning cluster

# In[475]:


resp_SF = resp_ad_cell_stim[:, 30:]
resp_SF.shape
resp_SF = norm_maxmin_across_stim(resp_SF)
model_SF = cluster_stim_resp(2, resp_SF, 5, 'grating and noise') # larger distance threshold


# In[476]:


plot_cluster_tuning(model_SF, resp_SF, 'SF', cluster_sortby='tuning peak', sharexy=False) # , ncol=5, figsize=(18, 30)


# ### all stim tuning cluster

# In[478]:


resp_all = resp_ad_cell_stim
resp_all = norm_maxmin_across_stim(resp_all)
model_all = cluster_stim_resp(3, resp_all, 5, 'all 50 stims') # larger distance threshold to avoid over clustering


# In[479]:


plot_cluster_tuning(model_all, resp_all, 'all', cluster_sortby='tuning peak', sharexy=False, ncol=3, figsize=(18, 30))
# # tight_layout cannot make axes height small enough to accommodate all axes decorations.
# # TODO: overlay tuning curve of grat, noise, nat, from same starting point


# In[274]:


# # imshow resp_ad_cell_stim after sorting cells by cluster

# # z score resp_ad_cell_stim
# resp_ad_cell_stim_sorted = (resp_ad_cell_stim - np.nanmean(resp_ad_cell_stim, axis=1, keepdims=True)) / np.nanstd(resp_ad_cell_stim, axis=1, keepdims=True)
# resp_ad_cell_stim_sorted = resp_ad_cell_stim_sorted[model_all.labels_, :]
# plt.figure(figsize=(10, 10))
# plt.imshow(resp_ad_cell_stim_sorted, aspect='auto', cmap='viridis')
# plt.colorbar()
# plt.xlabel('stim id')
# plt.ylabel('cell id')
# # plt.title('resp_ad_cell_stim_sorted')
# plt.tight_layout()


# ## clustering by adaptation

# In[480]:


resp_ad_by_stim[0].shape, len(resp_ad_by_stim) # ncell x ntrial_rep x nstim


def trial_rep_select(arr, n, seed=42):
    '''
    select 20 trials out of trial_rep
    arr: array to select from. shape [ncell x trial_rep]
    n: number of samples to select
    '''
    ntrial = arr.shape[1]
    np.random.seed(seed) # make sure R1 and R2 select the same trials!!
    idx = np.random.choice(ntrial, n, replace=True)
    return arr[:,idx]

    
def calc_adp_pop_boot_by_stim(resp_ad_by_stim, nboot=1000):
    '''calculate population adp by stim, with bootstrap sampling from trial reps'''
    
    adp_stim_boot = np.zeros((nboot, len(resp_ad_by_stim), resp_ad_by_stim[0].shape[0]))
    for iboot in np.arange(nboot):
        R1_agg_stim = np.array([trial_rep_select(irep, 20, seed=iboot).sum(axis=1) for irep in resp_ad_by_stim]) # sum over 20 trial reps
        R2_agg_stim = np.array([trial_rep_select(irep, 20, seed=iboot).sum(axis=1) for irep in resp_tg_by_stim])
        adp_stim = (R2_agg_stim - R1_agg_stim) / (R1_agg_stim + 1e-7)

        drop_perc = 0.5
        if cell_filter_name == 'vis':
            drop_perc = 3
        adp_stim[adp_stim < np.nanpercentile(adp_stim, drop_perc)] = np.nan
        adp_stim[adp_stim > np.nanpercentile(adp_stim, 100-drop_perc)] = np.nan # drop extreme values of adp_stim (top and bottom x%)
        adp_stim_boot[iboot] = adp_stim
    return adp_stim_boot

adp_stim_boot = calc_adp_pop_boot_by_stim(resp_ad_by_stim, nboot=1000)
adp_stim_boot.shape # nboot x nstim x ncell


# In[481]:


adp_cell_stim = adp_stim_boot.mean(axis=0).T
adp_cell_stim.shape # ncell x nstim


# ### grating adp cluster

# In[482]:


adp_grat = adp_cell_stim[:, 30:40]
# adp_grat = norm_adp_across_stim(adp_grat)
adp_grat = np.nan_to_num(adp_grat) # adp_grat set nan to 0 for clustering
model_grat_adp = cluster_stim_resp(0.5, adp_grat, 5, 'gratings adp')
plot_cluster_tuning(model_grat_adp, adp_grat, 'gratings', sharexy=False)


# In[ ]:


# TODO: take cluster labels from model_grat and plot adp_grat for each cluster. also try yuansi regression adp
# expectation: from current plots, looks like grating cluster will not have similar adp tuning curves

