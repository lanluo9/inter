#!/usr/bin/env python
# coding: utf-8

# In[6]:


os.chdir('d:\\repo\\inter_data\\inter\\')


# In[8]:


import numpy as np
import pandas as pd

import seaborn as sns
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

import os
from os import listdir
from os.path import join, isdir, isfile
import pickle

from src import adp

get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')


# In[9]:


mousedir = r'D:\repo\inter_data\mix50'.replace('\\', '/')

#Create a dictionary with keys being recording folder name, and value -- array with first element being stimulus identity array, and second element -- trace by trial matrix; 
recordings = {f : [adp.load_trace_trial_data(join(mousedir, f), vis_filter=True)] for f in listdir(mousedir) if isdir(join(mousedir, f))}


# In[10]:


# Calculate the minimum number of trials for each trial id across recordings and store in a dictionary with trial id as the key, and smallest number of trials as a value
min_trial_num = {}
for recording in recordings.keys():
    si, tbt = recordings[recording][0]
    names, counts = np.unique(si, return_counts=True)
    for i, n in enumerate(names):
        if n not in min_trial_num.keys() or counts[i] < min_trial_num[n]:
            min_trial_num[n] = counts[i]

#out_tbt is the final merged array
out_tbt = np.array([])
trialids = min_trial_num.keys()
# print(trialids)

for recording in recordings.keys():
    #for each recording, re-create trace-by-trial matrix with trial ids sorted the same way across recordings
    out_tbt_by_recording = np.array([])
    for trial_id in trialids:
        min_num_trials = min_trial_num[trial_id]
        si, tbt = recordings[recording][0]
        #select the trials where the current trial_id is used and use that to index
        curr_trialid_locs = np.where(si == trial_id)[:min_num_trials - 1][0]
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

out_tbt.shape


# In[14]:


# min_trial_num

stim_id_merged = np.array([])
for trial_id in trialids:
    min_num_trials = min_trial_num[trial_id]
    for i in range(min_num_trials):
        stim_id_merged = np.append(stim_id_merged, trial_id)

# count stim_id_merged elements
np.unique(stim_id_merged, return_counts=True)

