#!/usr/bin/env python
# coding: utf-8

# # import

# In[1]:


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


# In[ ]:


import scipy.io.loadmat


# In[2]:


local_flag = False
if local_flag:
    repo_dir = r'D:\repo\inter_data\inter'.replace("\\", "/") # under env dimred
else:
    repo_dir = r'C:\Users\ll357\Documents\inter'.replace("\\", "/")
os.chdir(repo_dir)
from src import adp


# # load grat1

# In[3]:


data_dir_grat = r'\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images\mat_inter\V1_i1372_220714_cellpose'
stim_id, trace_by_trial = adp.load_trace_trial_data(data_dir_grat, vis_filter=False)
stim_id = stim_id[:trace_by_trial.shape[1]] # truncate to uncorrupted trials
trace_avg_cell, trace_cell_sem, trace_stim_avg = adp.calc_trace_stim(trace_by_trial, stim_id)


# In[4]:


fig, ax = plt.subplots(figsize=(10, 5))
plt.plot(trace_avg_cell[:50], alpha=0.7, linewidth=3)
ax.axvspan(10, 12, alpha=0.2, color='red')
ax.axvspan(13, 15, alpha=0.2, color='gray')
ax.axvspan(20, 22, alpha=0.2, color='red')
plt.xlim(0, 50);


# # adp stability
# estimate how many trials to group together to get a stable adaptation value  
# 1. trace_by_trial shape: ncell x ntrial x nframe  
# 2. for each cell x trial, get resp_ad and resp_tg  
# 3. for whole population, loop thru trial numbers to group together to get a stable adaptation value  
# 3.1 measure adp value stability: plot group_size vs adp_mean & adp_std  
# 4. do this for each cell  
# 4.1 distribution of stable_group_size among cells  

# In[5]:


# get resp_ad and resp_tg [cell x trial]
base1 = trace_by_trial[:,:,0:2+1].mean(axis=2) # avg over time window frames
resp_ad = trace_by_trial[:,:,10:12+1].mean(axis=2)
resp_ad = resp_ad - base1

base2 = trace_by_trial[:,:,13:15+1].mean(axis=2)
resp_tg = trace_by_trial[:,:,20:22+1].mean(axis=2)
resp_tg = resp_tg - base2


# In[6]:


plt.hist(resp_ad.ravel(), bins=100, alpha=0.3, label='AD');
plt.hist(resp_tg.ravel(), bins=100, alpha=0.3, label='TG');
plt.xlim(-0.2, 0.6);


# In[10]:


# # cap resp to thres_perc percentile of trials

# thres_perc = 1
# low, high = np.percentile(resp_ad.flatten(), [thres_perc, 100-thres_perc])
# resp_ad = np.clip(resp_ad, low, high)
# plt.hist(resp_ad.ravel(), bins=100, alpha=0.3, label='AD');

# low, high = np.percentile(resp_tg.flatten(), [thres_perc, 100-thres_perc])
# resp_tg = np.clip(resp_tg, low, high)
# plt.hist(resp_tg.ravel(), bins=100, alpha=0.3, label='TG');


# ## pop adp stability
# population of 56 cells need 10-20 trials to get a stable adaptation value  
# regardless of using mean or median to aggregate  
# but using median gets abs larger adp value than mean

# In[7]:


mean_or_median = np.median # for population, median is better than mean

resp_ad_pop = mean_or_median(resp_ad, axis=0) # population response as avg or median of all cells for each trial
resp_tg_pop = mean_or_median(resp_tg, axis=0)

adp_norm_diff = (resp_tg_pop - resp_ad_pop) / (resp_ad_pop + 1e-7) # difference over baseline
adp_IOU = (resp_tg_pop - resp_ad_pop) / (resp_tg_pop + resp_ad_pop + 1e-7) # difference over sum
adp_norm_diff.shape, mean_or_median(adp_norm_diff), mean_or_median(adp_IOU) # group_size = 1, use all cell resp mean or median


# In[8]:


# shuffle resp_ad_pop and resp_tg_pop the same way
np.random.seed(42)
idx = np.random.permutation(resp_ad_pop.shape[0])
resp_ad_pop_shuf = resp_ad_pop[idx]
resp_tg_pop_shuf = resp_tg_pop[idx]

# for whole population, loop thru trial numbers to group together to get a stable adaptation value  
adp_agg = []
adp_std = []
for group_size in np.arange(1, trace_by_trial.shape[1] // 2):
    ngroup = trace_by_trial.shape[1] // group_size

    resp_ad_cut = resp_ad_pop_shuf[:group_size*ngroup].reshape(ngroup, group_size) # reshape to ngroup x group_size
    resp_tg_cut = resp_tg_pop_shuf[:group_size*ngroup].reshape(ngroup, group_size)
    resp_ad_group = mean_or_median(resp_ad_cut, axis=1) # aggregate resp within group of trials
    resp_tg_group = mean_or_median(resp_tg_cut, axis=1)

    adp_group = (resp_tg_group - resp_ad_group) / (resp_ad_group + 1e-7) # calc adp with `avg resp within group`
    adp_group_agg = mean_or_median(adp_group) # agg adp across group
    adp_group_std = np.std(adp_group) # std of adp across group
    adp_agg.append(adp_group_agg) # group_size variable, all cell resp agg
    adp_std.append(adp_group_std)

adp_agg = np.array(adp_agg)
adp_std = np.array(adp_std)
adp_sem = adp_std / np.sqrt(adp_std.shape[0])


# In[9]:


plt.figure(figsize=(12, 10))
plt.subplot(2,1,1)
plt.plot(adp_agg, 'o-', color='b', linewidth=3, alpha=0.3, label='adp agg across groups')
# plt.axhline(-0.32, color='g')
plt.xticks(np.arange(0,110,10))
plt.ylim(plt.ylim()[0], 0)
plt.legend(frameon=False);

plt.subplot(2,1,2)
plt.plot(adp_std, 'o-', color='r', linewidth=3, alpha=0.5, label='adp std across groups')
plt.xticks(np.arange(0,110,10))
plt.xlabel('group_size - 1')
# plt.ylim(0, 0.5)
plt.ylabel('adp value')
plt.legend(frameon=False);


# In[10]:


errbar = adp_std # adp_sem 
x = np.arange(1, adp_agg.shape[0]+1)

plt.style.use('seaborn-white')

fig = plt.figure(figsize=(12, 6))
plt.plot(x, adp_agg, '-', color='b', linewidth=3, alpha=0.3)
plt.fill_between(x, 
                 adp_agg + errbar, adp_agg - errbar, 
                 color='r', alpha=0.2)
plt.axhline(0, color='k', linewidth=1, linestyle='--')

plt.xticks(np.arange(0,110,10))
plt.xlabel('number of trials grouped together', fontsize=18)
plt.ylabel('adaptation index of population across groups', fontsize=18); #  across groups of different sizes


# ## single cell adp stability
# single cells need 20-30 trials to get a stable adaptation value  

# In[11]:


mean_or_median = np.median # for single cell, median is better than mean too

ncell = trace_by_trial.shape[0]
ngroup_min = 2
nsize = len(np.arange(1, trace_by_trial.shape[1] // ngroup_min)) # how many group_size to test
adp_cell_agg = np.zeros((ncell, nsize))
adp_cell_std = np.zeros((ncell, nsize))
adp_cell_sem = np.zeros((ncell, nsize))

for icell in np.arange(ncell):
    resp_ad_cell = resp_ad[icell, :]
    resp_tg_cell = resp_tg[icell, :]

    # shuffle resp_ad_pop and resp_tg_pop the same way
    np.random.seed(42)
    idx = np.random.permutation(resp_ad_pop.shape[0])
    resp_ad_cell_shuf = resp_ad_cell[idx]
    resp_tg_cell_shuf = resp_tg_cell[idx]

    for isize, group_size in enumerate(np.arange(1, trace_by_trial.shape[1] // ngroup_min)):
        ngroup = trace_by_trial.shape[1] // group_size

        resp_ad_cut = resp_ad_cell_shuf[:group_size*ngroup].reshape(ngroup, group_size)
        resp_tg_cut = resp_tg_cell_shuf[:group_size*ngroup].reshape(ngroup, group_size)
        resp_ad_group = np.sum(resp_ad_cut, axis=1)
        resp_tg_group = np.sum(resp_tg_cut, axis=1)

        adp_group = (resp_tg_group - resp_ad_group) / (resp_ad_group + 1e-7)
        adp_group_agg = mean_or_median(adp_group)
        adp_group_std = np.std(adp_group)
        adp_group_sem = adp_group_std / np.sqrt(adp_group.shape[0])
        adp_cell_agg[icell, isize] = adp_group_agg
        adp_cell_std[icell, isize] = adp_group_std
        adp_cell_sem[icell, isize] = adp_group_sem


# In[12]:


adp_cell_std_agg = mean_or_median(adp_cell_std, axis=0)

plt.figure(figsize=(12, 5))
plt.plot(adp_cell_std_agg, 'o-', color='r', linewidth=3, alpha=0.7, label='adp std across groups, agg over cells')
for icell in np.arange(ncell):
    plt.plot(adp_cell_std[icell, :], '-', color='r', linewidth=3, alpha=0.05, label='adp std across groups')

plt.axhline(1, alpha=0.4, color='gray')
plt.xticks(np.arange(0,110,10))
plt.xlabel('group_size - 1')

# plt.ylim(1e-3, plt.ylim()[1]) # set y axis lower limit
# plt.ylim(1e-2, 1e2)
plt.yscale('log')
plt.ylabel('adp std');


# In[13]:


adp_cell_agg_agg = mean_or_median(adp_cell_agg, axis=0)

plt.figure(figsize=(12, 5))
plt.plot(adp_cell_agg_agg, 'o-', color='b', linewidth=3, alpha=0.7, label='adp agg across groups, agg over cells')
# for icell in np.arange(ncell):
#     plt.plot(adp_cell_agg[icell, :], '-', color='b', linewidth=3, alpha=0.05, label='adp agg across groups')

plt.xticks(np.arange(0,110,10))
plt.xlabel('group_size - 1')

# plt.yscale('symlog') # plot negative values on -log scale, plot ~0 values linearly: https://stackoverflow.com/questions/3305865/what-is-the-difference-between-log-and-symlog
# plt.ylim(-0.5,1)
plt.ylabel('adp mean');


# In[14]:


plt.plot(x, adp_agg, '-', color='b', linewidth=3, alpha=0.3)
plt.fill_between(x, 
                 adp_agg + errbar, adp_agg - errbar, 
                 color='r', alpha=0.2)
plt.axhline(0, color='k', linewidth=1, linestyle='--')

plt.xticks(np.arange(0,110,10))
plt.xlabel('number of trials grouped together', fontsize=18)
plt.ylabel('adaptation index of population', fontsize=18);


# In[20]:


adp_cell_sem_agg = mean_or_median(adp_cell_sem, axis=0) # adp error between groups -> agg across cells
errbar = adp_cell_sem_agg
x = np.arange(1, adp_agg.shape[0]+1)
xlim = 50

plt.rcParams.update(plt.rcParamsDefault)

fig = plt.figure(figsize=(6, 4))
plt.plot(x[:xlim], adp_cell_agg_agg[:xlim], '-', color='blue', linewidth=3, alpha=0.7, label='single neurons')
plt.fill_between(x[:xlim],  
                adp_cell_agg_agg[:xlim] + errbar[:xlim], adp_cell_agg_agg[:xlim] - errbar[:xlim], color='b', alpha=0.2)

errbar = adp_sem 
# x = np.arange(1, adp_agg.shape[0]+1)
plt.plot(x[:xlim], adp_agg[:xlim], '-', color='red', linewidth=3, alpha=0.3, label='population')
plt.fill_between(x[:xlim], 
                 adp_agg[:xlim] + errbar[:xlim], adp_agg[:xlim] - errbar[:xlim], 
                 color='r', alpha=0.2)

plt.axhline(0, color='k', linewidth=1, linestyle='--')

plt.xticks(np.arange(0, xlim+1, 10), fontsize=16)
plt.yticks(fontsize=16)
plt.xlabel('Number of trials grouped together', fontsize=18)
plt.ylabel('Adaptation index', fontsize=18)
plt.legend(fontsize=18, loc='lower right', frameon=False);

# top and right spines off
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.show()
plt.draw()
# figpath = r'C:\Users\ll357\Documents\inter\results\poster SFN 2022/'.replace('\\', '/')
# fig.savefig(os.path.join(figpath, 'adp_stability_group_trials.pdf'), format='pdf')

fig4d_dict = {'xlim': xlim, 'x': x, 'adp_cell_agg_agg': adp_cell_agg_agg, 'adp_cell_sem_agg': adp_cell_sem_agg, 'adp_agg': adp_agg, 'adp_sem': adp_sem}
with open(os.path.join('fig4d_dict.pkl'), 'wb') as f:
    pickle.dump(fig4d_dict, f)


# # load mix50
# mix50 stim: Z:\All_Staff\home\lan\Mwork\mix50 - bunnytop high lum contrast mix grating and noise\Image

# In[22]:


mousedir = r'D:\repo\inter_data\mix50'.replace('\\', '/')
# imouse = 'i1373'
# date = '220909'
# area = 'V1'
# if local_flag:
#     data_dir = r'D:\repo\inter_data\mix50\V1_i1373_220909_cellpose'
# else:
#     data_dir = r'\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images\mat_inter/' + area + '_' + imouse + '_' + date + '_cellpose/'

# Create a dictionary with keys being recording folder name, and value -- array with first element being stimulus identity array, and second element -- trace by trial matrix; 
recordings = {f : [adp.load_trace_trial_data(os.path.join(mousedir, f), vis_filter=True)] for f in os.listdir(mousedir) if os.path.isdir(os.path.join(mousedir, f))}

# Calculate the minimum number of trials for each trial id across recordings and store in a dictionary with trial id as the key, and smallest number of trials as a value
min_trial_num = {}
for recording in recordings.keys():
    si, tbt = recordings[recording][0]
    names, counts = np.unique(si, return_counts=True)
    for i, n in enumerate(names):
        if n not in min_trial_num.keys() or counts[i] < min_trial_num[n]:
            min_trial_num[n] = counts[i]


# In[152]:


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


# In[153]:


stim_id_merged = np.array([])
for trial_id in trialids:
    min_num_trials = min_trial_num[trial_id]
    for i in range(min_num_trials):
        stim_id_merged = np.append(stim_id_merged, int(trial_id))
stim_id_merged = stim_id_merged.astype(int)

# count stim_id_merged elements
np.unique(stim_id_merged, return_counts=True)
stim_id_merged.shape


# In[154]:


trace_by_trial = out_tbt
stim_id = stim_id_merged
# stim_id, trace_by_trial = adp.load_trace_trial_data(data_dir, vis_filter=True)
# print(trace_by_trial.shape)
trace_avg_cell, trace_cell_sem, trace_stim_avg = adp.calc_trace_stim(trace_by_trial, stim_id)

fig, ax = plt.subplots(figsize=(10, 5))
plt.plot(trace_avg_cell[:50], alpha=0.7, linewidth=3)
ax.axvspan(9,11, alpha=0.2, color='red')
ax.axvspan(13,14, alpha=0.2, color='gray')
ax.axvspan(20,22, alpha=0.2, color='red')
plt.xlim(0, 50);


# # trace by stim / type

# In[155]:


# plot trace for each stim
plt.figure(figsize=(20, 7))
plt.plot(trace_avg_cell, alpha=0.8, linewidth=7, color='cyan')
for i in np.unique(stim_id)-1:
    plt.plot(trace_stim_avg[i], alpha=0.2)


# In[156]:


# plot trace for stim 1-30 (natural) vs stim 31-40 (grat) vs stim 41-50 (noise)

stim_type_dict = {'natural': [np.arange(1, 30+1)], 
                  'grat': [np.arange(31, 40+1)], 
                  'noise': [np.arange(41, 50+1)]}
for key in stim_type_dict.keys():
    print(key, np.min(stim_type_dict[key]), '-', np.max(stim_type_dict[key]))
    
trace_type_avg = []
trace_type_std = []
trace_type_sem = []

for key in stim_type_dict.keys():
    
    trace_itype_avg = np.mean(trace_by_trial[:, np.where((np.isin(stim_id, stim_type_dict[key])))[0]], axis=1) # ncell x nframe
    trace_itype_avg = np.mean(trace_itype_avg, axis=0) # nframe
    trace_itype_std = np.std(trace_by_trial[:, np.where((np.isin(stim_id, stim_type_dict[key])))[0]], axis=1)
    trace_itype_std = np.mean(trace_itype_std, axis=0)
    trace_itype_sem = trace_itype_std / np.sqrt(len(np.where((np.isin(stim_id, stim_type_dict[key])))[0]))

    trace_type_avg.append(trace_itype_avg)
    trace_type_std.append(trace_itype_std)
    trace_type_sem.append(trace_itype_sem)

len(trace_type_avg), trace_type_avg[0].shape


# In[157]:


# plot trace for each stim type
plt.figure(figsize=(10, 6))

for i in np.arange(len(stim_type_dict.keys())):
    plt.plot(trace_type_avg[i], alpha=0.5, linewidth=3, label=list(stim_type_dict.keys())[i])
    errbar = trace_type_sem[i]
    plt.fill_between(np.arange(len(trace_type_avg[i])), 
                 trace_type_avg[i] + errbar, trace_type_avg[i] - errbar, 
                 alpha=0.1, label='')

plt.xlim(0, 50)
plt.legend(frameon=False);


# In[ ]:


# poster SFN fig3a

plt.figure(figsize=(6, 5))

for i in np.arange(len(stim_type_dict.keys())-1): # only plot natural and grat
    i = 1-i # flip order to: first grat, then natural
    plt.plot(trace_type_avg[i], alpha=0.5, linewidth=3, label=list(stim_type_dict.keys())[i])
    errbar = trace_type_sem[i]
    plt.fill_between(np.arange(len(trace_type_avg[i])), 
                 trace_type_avg[i] + errbar, trace_type_avg[i] - errbar, 
                 alpha=0.1, label='')

plt.xlim(0, 70)
plt.legend(frameon=False);


# ## poster SFN fig3a

# In[ ]:


# plt default rcparams
plt.rcParams.update(plt.rcParamsDefault)

plt.figure(figsize=(6, 5))

plt.figure(figsize=(6, 5))

for i in np.arange(len(stim_type_dict.keys())-1): # only plot natural and grat
    i = 1-i # flip order to: first grat, then natural
    plt.plot(trace_type_avg[i], alpha=0.5, linewidth=3, label=list(stim_type_dict.keys())[i])
    errbar = trace_type_sem[i]
    plt.fill_between(np.arange(len(trace_type_avg[i])), 
                 trace_type_avg[i] + errbar, trace_type_avg[i] - errbar, 
                 alpha=0.1, label='')

# change xticks to time in sec
frame_rate = 30 # frame per second
nlabel = 12
x = np.arange(0, trace_type_avg[i].shape[0], nlabel)
time_sec = np.arange(0, trace_type_avg[i].shape[0], nlabel) / frame_rate
labels = [str(t) for t in time_sec]
plt.xticks(x, labels, fontsize=16);
plt.yticks(fontsize=16)
plt.xlim([0,70])

plt.xlabel('Time (s)', fontsize=18)
plt.ylabel('dF/F', fontsize=18)
plt.grid(False)
plt.legend(frameon=False, fontsize=16);

ax = plt.gca() # spine off for right and top
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

plt.tight_layout()
# plt.savefig('trace_grand_avg_compare_stim_type.pdf')


# # resp ad vs tg

# In[158]:


# get resp_ad and resp_tg [cell x trial]
base1 = trace_by_trial[:,:,0:2+1].mean(axis=2) # avg over time window frames
resp_ad = trace_by_trial[:,:,9:11+1].mean(axis=2)
resp_ad = resp_ad - base1

base2 = trace_by_trial[:,:,13:14+1].mean(axis=2)
resp_tg = trace_by_trial[:,:,20:22+1].mean(axis=2)
resp_tg = resp_tg - base2

resp_ad.shape, resp_tg.shape


# # feature matrix for decoding branch
# did not edit to accomodate mix50 merged data yet

# In[384]:


feature_mat = np.concatenate((resp_ad.T, resp_tg.T), axis=0)
feature_mat.shape # ntrial x ncell. doubled ntrial bc we now count ad and tg as separate trials

feature_df = pd.DataFrame(feature_mat)
feature_df.columns = [f'{date}_cell_' + str(i) for i in np.arange(feature_mat.shape[1])] # name columns

image_id = stim_id + stim_id # pass thru adapter stim id and target stim id. they are identical
feature_df['image_index'] = image_id

feature_df['repeat_number'] = [0] * len(stim_id) + [1] * len(stim_id) # first show of an image is not a repeat, repeat_number = 0
feature_df['is_change'] = feature_df['repeat_number'].astype('bool')
feature_df['is_change'] = ~feature_df['is_change'] # repeat_number = 0 means is_change = True (R1)

feature_df


# In[386]:


feature_df.groupby(['repeat_number', 'image_index']).count()#.get_group((0, 1))


# In[10]:


plt.hist(resp_ad.ravel(), bins=100, alpha=0.3, label='AD');
plt.hist(resp_tg.ravel(), bins=100, alpha=0.3, label='TG');
plt.xscale('symlog')
plt.yscale('symlog')


# In[11]:


# # set extreme percentile resp to nan

# thres_perc_low = 0.5
# thres_perc_high = 0.5
# low, high = np.percentile(resp_ad.flatten(), [thres_perc_low, 100-thres_perc_high])
# # resp_ad = np.clip(resp_ad, low, high)
# resp_ad[resp_ad < low] = np.nan
# resp_ad[resp_ad > high] = np.nan
# plt.hist(resp_ad.ravel(), bins=100, alpha=0.3, label='AD');

# low, high = np.percentile(resp_tg.flatten(), [thres_perc, 100-thres_perc])
# # resp_tg = np.clip(resp_tg, low, high)
# resp_tg[resp_tg < low] = np.nan
# resp_tg[resp_tg > high] = np.nan
# plt.hist(resp_tg.ravel(), bins=100, alpha=0.3, label='TG')
# plt.legend();


# ## adp or resp corr btw stim

# In[432]:


R1_R2_avg = feature_df.groupby(['image_index', 'repeat_number']).mean().reset_index() # avg R1 & R2 for each img, each neuron
non_neuron_columns = ['image_index', 'repeat_number', 'is_change'] # in our data, cell id col name are also strings
R1_avg = R1_R2_avg[R1_R2_avg.repeat_number == 0].drop(non_neuron_columns, axis=1).reset_index(drop=True).to_numpy()
R2_avg = R1_R2_avg[R1_R2_avg.repeat_number == 1].drop(non_neuron_columns, axis=1).reset_index(drop=True).to_numpy()
adp_img = (R2_avg - R1_avg) / (R1_avg + 1e-7) # adp for each img (trial already grouped), each neuron
nstim = feature_df.image_index.unique().shape[0] # use trials_to_keep to get nstim, bc in full data exists img_id=8 where repeat_number=nan

def plot_sorted_cells(adp_img, title_str):
    # for each row in adp_img, mask cells with extreme values by percentile
    adp_img_masked = np.ma.masked_where(adp_img > np.percentile(adp_img, 99), adp_img)
    adp_img_masked = np.ma.masked_where(adp_img_masked < np.percentile(adp_img_masked, 1), adp_img_masked)
    print(np.amax(adp_img).round(2), np.amin(adp_img).round(2))
    print(np.amax(adp_img_masked).round(2), np.amin(adp_img_masked).round(2))

    adp_sorted_idx = np.argsort(adp_img_masked[0, :]) # cell sorted by adp to img 0
    adp_sorted_idx = adp_sorted_idx[::-1] # argsort from large to small, so reverse the order

    # make subplots for each img pair: img 0 vs img 1, img 0 vs img 2, etc
    fig, axs = plt.subplots(11, 5, figsize=(20, 30), sharex=True, sharey=True)
    plt.suptitle(title_str)
    adp_stim0 = adp_img_masked[0, adp_sorted_idx].round(2)
    corr_to_stim0 = np.ones((nstim, 1))

    for istim in np.arange(nstim):
        adp_stim = adp_img_masked[istim, adp_sorted_idx]
        corr = np.corrcoef(adp_stim0, adp_stim)[0, 1]
        color_str = 'k' if corr < 0.1 else 'r'
        corr_to_stim0[istim] = corr
        axs[istim // 5, istim % 5].text(0.5, 0.9, 'corr: ' + str(corr.round(2)), horizontalalignment='center', verticalalignment='center', transform=axs[istim // 5, istim % 5].transAxes, color=color_str, fontsize=12)
        
        axs[istim // 5, istim % 5].plot(adp_stim0, label=str(0), alpha=0.5, color='g', linewidth=3)
        axs[istim // 5, istim % 5].plot(adp_stim, label=str(istim), alpha=0.5)
        axs[istim // 5, istim % 5].set_title('img 0 & ' + str(istim))
    # TODO: sort cells by adp avg over all trials, then plot adp to each img

    adp_stim_other = np.ma.median(adp_img_masked[1:, adp_sorted_idx], axis=0)
    # print(np.round(adp_stim_other, 2))
    corr = np.corrcoef(adp_stim0, adp_stim_other)[0, 1]
    color_str = 'k' if corr < 0.1 else 'r'
    axs[10, 4].text(0.5, 0.9, 'corr: ' + str(corr.round(2)), horizontalalignment='center', verticalalignment='center', transform=axs[10, 4].transAxes, color=color_str, fontsize=14)
    axs[10, 4].plot(adp_stim0, label=str(0), alpha=0.5, color='g', linewidth=3)
    axs[10, 4].plot(adp_stim_other, label='other', alpha=0.5, color='orange')
    axs[10, 4].set_title('img 0 & median of other stims') # final subplot: median resp to all stim except 0
    # plt.tight_layout()

    return corr_to_stim0

corr_to_stim0_adp = plot_sorted_cells(adp_img, 'adp to each img')


# In[433]:


corr_to_stim0_R1 = plot_sorted_cells(R1_avg, 'R1 to each img')


# In[434]:


corr_to_stim0_R2 = plot_sorted_cells(R2_avg, 'R2 to each img')


# In[437]:


# migrated from https://github.com/lanluo9/dim_red_decoding/blob/master/code/dim_red.ipynb

# question: 
# are response patterns more separable post adaptation?
# do sparseness increase after adaptation?

# method:
# sort cells by response to image #1    
# plot response pattern to image #2  
# do this for: pre vs post adaptation, diff image pairs  
# compare: pearson correlation of sorted cell response vector

# result:
# no obvious decrease in correlation after adaptation
# TODO: check diff img pair combinations, check diff istim vs other, avg over all img pairs

plt.plot(corr_to_stim0_adp, label='adp')
plt.plot(corr_to_stim0_R1, label='R1')
plt.plot(corr_to_stim0_R2, label='R2')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False)


# # adp
# ## adp by stim type: agg over all cell all trial
# 
# previous adp def  
# 
# adp_norm_diff = (resp_tg_pop - resp_ad_pop) / (resp_ad_pop + 1e-7) # difference over baseline  
# adp_IOU = (resp_tg_pop - resp_ad_pop) / (resp_tg_pop + resp_ad_pop + 1e-7) # difference over sum  
# adp_norm_diff.shape, mean_or_median(adp_norm_diff), mean_or_median(adp_IOU) # group_size = 1, use all cell resp mean or median

# In[159]:


agg_fun = np.nansum

resp_ad_pop = agg_fun(resp_ad, axis=0) # population response as agg of all cells for each trial
resp_tg_pop = agg_fun(resp_tg, axis=0)

df = pd.DataFrame()
df['stim_id'] = stim_id
df['stim_type'] = np.where(np.isin(stim_id, stim_type_dict['natural']), 'natural',
                  np.where(np.isin(stim_id, stim_type_dict['grat']), 'grat',
                  np.where(np.isin(stim_id, stim_type_dict['noise']), 'noise', 'other')))

# print(len(df))
# print(resp_ad_pop.shape)
df['resp_ad_pop'] = resp_ad_pop
df['resp_tg_pop'] = resp_tg_pop
# df = df.sort_values(by=['stim_id'])
# df

df_pop = df.groupby('stim_type').sum() # agg pop resp of all trials for each stim type
df_pop['adp_pop_type'] = (df_pop['resp_tg_pop'] - df_pop['resp_ad_pop']) / (df_pop['resp_ad_pop'] + 1e-7)
df_pop


# ## vis, img driven, pref filter
# merged mix50 data only support 'vis' filter for now

# In[243]:


img_driven_merge = np.array([])

for subdir in os.listdir(r'D:\repo\inter_data\mix50'):
    data_dir = os.path.join(r'D:\repo\inter_data\mix50', subdir)
    print(data_dir)

    with open(data_dir + "/vis_driven.pickle", "rb") as handle:
        vis = pickle.load(handle)
    vis_driven = vis["vis_driven"]
    vis_driven = [v[0] for v in vis_driven]

    img_driven = vis['img_driven'][vis_driven, :]
    print(img_driven.shape)
    img_driven_merge = img_driven if img_driven_merge.size == 0 else np.concatenate((img_driven_merge, img_driven))

img_driven = img_driven_merge
img_driven.shape


# In[244]:


# associate w spatial freq after bootstrap. calc pop resp w img_driven filter

# # data_dir = r'\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images\mat_inter\V1_i1373_220909_cellpose'
# data_dir = r'D:\repo\inter_data\mix50\V1_i1373_220909_cellpose'
# with open(data_dir + "/vis_driven.pickle", "rb") as handle:
#     vis = pickle.load(handle)
# # vis.keys()

# vis_driven = vis["vis_driven"]
# vis_driven = [v[0] for v in vis_driven]
# img_driven = vis['img_driven'][vis_driven, :] # filter out non-vis driven cells. match shape of resp_ad
# # img_driven.shape

# # do img driven cells add up to vis driven cells? no.
# np.sum(img_driven, axis=1) # nimg that drive each cell
# ncell_img_driven = np.sum(np.sum(img_driven, axis=1) > 0) # ncell that are img driven
# ncell_img_driven < np.sum(vis_driven) # some vis driven cells are not img driven! 
# # TODO: adjust vis_driven criteria: lindsey uses union of img driven cells

cell_filter_name = 'vis'
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


# In[245]:


plt.matshow(cell_filter.T)
plt.xticks([]) # no axis ticks
plt.yticks([])
plt.xlabel('cells')
plt.ylabel('stims');

print(cell_filter.shape) # ncell x nstim
print(f'{np.sum(cell_filter.sum(axis=0) == 0)} out of {cell_filter.shape[1]} stims have no cells passing filter')
print(f'{sum(cell_filter.sum(axis=1) == 0)} out of {cell_filter.shape[0]} cells have no stim passing filter')


# ## apply filter

# In[246]:


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
len(resp_ad_by_stim), [i.shape for i in resp_ad_by_stim[:5]]


# ## adp w bootstrap sampling

# ### boot adp by stim type

# In[247]:


## bootstrap 20 out of 30 trials
'''
for each stim id, get resp_ad and resp_tg [ncell x trial_rep x nstim]
bootstrap select 20 trials out of trial_rep
sum over 20 trials to get resp & adp [ncell x nstim]
'''

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
    
# sum over 20 trials to get resp & adp [nstim x ncell]

seed = 42
R1_agg_stim = np.array([trial_rep_select(i, 20, seed=seed).sum(axis=1) for i in resp_ad_by_stim]) # sum over 20 trial reps
R2_agg_stim = np.array([trial_rep_select(i, 20, seed=seed).sum(axis=1) for i in resp_tg_by_stim])

adp_stim = (R2_agg_stim - R1_agg_stim) / (R1_agg_stim + 1e-7)
drop_perc = 0.5
if cell_filter_name == 'vis':
    drop_perc = 5
adp_stim[adp_stim < np.nanpercentile(adp_stim, drop_perc)] = np.nan
adp_stim[adp_stim > np.nanpercentile(adp_stim, 100-drop_perc)] = np.nan # drop extreme values of adp_stim (top and bottom x%)
print(adp_stim.shape)

plt.hist(adp_stim[:30, :].ravel(), bins=20, alpha=0.3, label='natural', density=True);
plt.hist(adp_stim[30:40, :].ravel(), bins=20, alpha=0.3, label='grat', density=True);
plt.hist(adp_stim[40:, :].ravel(), bins=20, alpha=0.3, label='noise', density=True);

plt.axvline(np.nanmedian(adp_stim[:30, :].ravel()), color='C0', label='median natural')
plt.axvline(np.nanmedian(adp_stim[30:40, :].ravel()), color='C1', label='median grat')
plt.axvline(np.nanmedian(adp_stim[40:, :].ravel()), color='C2', label='median noise')
plt.legend(frameon=False);


# ### viz resp of single cell by stim type

# In[248]:


R1_agg_stim.shape # nstim x ncell

df_resp_nat = pd.DataFrame(R1_agg_stim[:30, :].ravel(), columns=['resp'])
df_resp_nat['stim_type'] = 'natural'
df_resp_grat = pd.DataFrame(R1_agg_stim[30:40, :].ravel(), columns=['resp'])
df_resp_grat['stim_type'] = 'grat'
df_resp_noise = pd.DataFrame(R1_agg_stim[40:, :].ravel(), columns=['resp'])
df_resp_noise['stim_type'] = 'noise'
df_resp_stim_type = pd.concat([df_resp_nat, df_resp_grat, df_resp_noise], axis=0)
df_resp_stim_type['resp_number'] = 'R1'

df_R2_stim_type = df_resp_stim_type.copy()
df_R2_stim_type['resp'] = R2_agg_stim[:, :].ravel()
df_R2_stim_type['resp_number'] = 'R2'
df_resp_stim_type = pd.concat([df_resp_stim_type, df_R2_stim_type], axis=0)

df_resp_stim_type = df_resp_stim_type.dropna().reset_index(drop=True)
df_resp_stim_type


# In[249]:


# # boxplot for df_adp_stim_type
# sns.boxplot(x='stim_type', y='adp', data=df_adp_stim_type, palette='Set2')
# plt.ylabel('adp');

# violinplot for df_adp_stim_type
sns.violinplot(x='stim_type', y='resp', data=df_resp_stim_type, 
               hue="resp_number", split=True,
               palette='Set3', 
               inner="quartile",
               )
plt.ylabel('aggregated R1 over trial group');
plt.ylim(-2, 8)
plt.legend(frameon=False, );

# # swarmplot for df_adp_stim_type
# plt.figure(figsize=(20, 4))
# sns.swarmplot(x='stim_type', y='adp', data=df_adp_stim_type, palette='Set2')
# plt.ylabel('adp');


# In[250]:


adp_stim.shape # nstim x ncell

df_adp_stim_nat = pd.DataFrame(adp_stim[:30, :].ravel(), columns=['adp'])
df_adp_stim_nat['stim_type'] = 'natural'
df_adp_stim_grat = pd.DataFrame(adp_stim[30:40, :].ravel(), columns=['adp'])
df_adp_stim_grat['stim_type'] = 'grating'
df_adp_stim_noise = pd.DataFrame(adp_stim[40:, :].ravel(), columns=['adp'])
df_adp_stim_noise['stim_type'] = 'noise'

df_adp_stim_type = pd.concat([df_adp_stim_nat, df_adp_stim_grat, df_adp_stim_noise])
df_adp_stim_type = df_adp_stim_type.dropna().reset_index(drop=True)
df_adp_stim_type['cell_filter_name'] = 'img_driven'
df_adp_stim_type


# In[45]:


# df_adp_stim_nat = pd.DataFrame(adp_stim[:30, :].ravel(), columns=['adp'])
# df_adp_stim_nat['stim_type'] = 'natural'
# df_adp_stim_grat = pd.DataFrame(adp_stim[30:40, :].ravel(), columns=['adp'])
# df_adp_stim_grat['stim_type'] = 'grating'
# df_adp_stim_noise = pd.DataFrame(adp_stim[40:, :].ravel(), columns=['adp'])
# df_adp_stim_noise['stim_type'] = 'noise'

# df_adp_stim_type_pref = pd.concat([df_adp_stim_nat, df_adp_stim_grat, df_adp_stim_noise])
# df_adp_stim_type_pref = df_adp_stim_type_pref.dropna().reset_index(drop=True)
# df_adp_stim_type_pref['cell_filter_name'] = 'pref_filter'
# df_adp_stim_type_pref


# In[46]:


# df_adp_stim_type = pd.concat([df_adp_stim_type, df_adp_stim_type_pref])
# df_adp_stim_type


# ### viz boot adp for single cell by stim type

# In[251]:


# # boxplot for df_adp_stim_type
# sns.boxplot(x='stim_type', y='adp', data=df_adp_stim_type, palette='Set2')
# plt.ylabel('adp');

# violinplot for df_adp_stim_type
sns.violinplot(x='stim_type', y='adp', data=df_adp_stim_type, 
            #    hue="cell_filter_name", split=True,
               palette='Set3', 
               inner="quartile",
               )
plt.ylabel('adaptation'); # calc from bootstraped single cell resp
# plt.ylim(-2, 2)
plt.legend(frameon=False, );
print('median natural: ', np.nanmedian(adp_stim[:30, :].ravel()))
print('median grat: ', np.nanmedian(adp_stim[30:40, :].ravel()))
print('median noise: ', np.nanmedian(adp_stim[40:, :].ravel()))

# # swarmplot for df_adp_stim_type
# plt.figure(figsize=(20, 4))
# sns.swarmplot(x='stim_type', y='adp', data=df_adp_stim_type, palette='Set2')
# plt.ylabel('adp');


# In[252]:


# resp matrix [nstim x ncell] w filter

plt.figure(figsize=(20,8))
plt.subplot(2,1,1)
plt.matshow(R1_agg_stim, fignum=False)
plt.xticks([]) # no axis ticks
plt.yticks([]);
plt.colorbar()

plt.subplot(2,1,2)
plt.matshow(adp_stim, fignum=False)
plt.xticks([]) # no axis ticks
plt.yticks([]);
plt.colorbar()
plt.tight_layout()


# ### boot adp for pop

# In[253]:


resp_ad_by_stim[0].shape, len(resp_ad_by_stim) # ncell x ntrial_rep x nstim

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


# ### boot adp pop vs spatial freq

# In[254]:


cpd_array = np.round(np.geomspace(0.03, 0.9, num=10), 2)
noise_nbyn = [2,3,5,17, 6,10,34, 15,51,85] # from do_image_preprocess.py -> 2x2 grid to 85x85 grid
spatial_freq_cpd = list(cpd_array)

adp_stim_boot_pop = np.nanmean(adp_stim_boot, axis=2) # [nboot x nstim], avg over cells
adp_stim_pop_mean = np.nanmean(adp_stim_boot_pop, axis=0) # [nstim], avg over boot
adp_stim_pop_err = np.nanstd(adp_stim_boot_pop, axis=0)
adp_stim_pop_CI95 = np.nanpercentile(adp_stim_boot_pop, [2.5, 97.5], axis=0) # [2 x nstim], 95% CI


# In[255]:


def plot_adp_by_SF(cell_filter, adp_stim_pop_mean, adp_stim_pop_CI95, adp_stim_pop_err):
    x = np.arange(len(adp_stim_pop_mean[30:]))
    y1 = adp_stim_pop_mean[30:]
    y2 = np.sum(cell_filter[:, 30:], axis=0)

    fig, ax1 = plt.subplots(figsize=(8,5))
    ax2 = ax1.twinx()

    ax1.plot(x, y1, 'g-')
    ax1.fill_between(x, y1-adp_stim_pop_err[30:], y1+adp_stim_pop_err[30:], alpha=0.3, color='g');
    ax1.fill_between(x, adp_stim_pop_CI95[0, 30:], adp_stim_pop_CI95[1, 30:], alpha=0.3, color='orange');
    ax1.axhline(0, color='k', linestyle='-', alpha=0.3)

    ax2.plot(x, y2, 'b-')
    ax2.axvline(len(x)//2 - 0.5, color='k', linestyle='-.', alpha=0.3)

    ax1.set_xlabel('spatial frequency (cpd)')
    ax1.set_xticks(x, list(spatial_freq_cpd) + list(noise_nbyn), rotation=45);
    ax1.set_ylabel('adaptation index', color='g')
    ax2.set_ylabel(f'filtered ncell: {cell_filter_name}', color='b')

    ax1.grid(True)
    ax2.grid(False) # grid off
    plt.tight_layout()
    plt.show()

plot_adp_by_SF(cell_filter, adp_stim_pop_mean, adp_stim_pop_CI95, adp_stim_pop_err)


# #### bin SF
# get more ncell preferring each SF

# In[256]:


stim_id_mod = stim_id.copy()

stim_id_mod = np.array(stim_id_mod) # convert to array, otherwise bool index will not work!

stim_id_mod[stim_id_mod==31] = 32 # pretend stim #31 is stim #32: merge low SF gratings
# stim_id_mod[stim_id_mod==39] = 38 
# stim_id_mod[stim_id_mod==40] = 38
stim_id_mod[stim_id_mod==40] = 39 # bin high SF gratings

stim_id_mod[stim_id_mod==41] = 42 # same for noise
# stim_id_mod[stim_id_mod==49] = 48
# stim_id_mod[stim_id_mod==50] = 48
stim_id_mod[stim_id_mod==50] = 49

np.unique(stim_id_mod), len(stim_id_mod) # stim id starts from 1!


# In[257]:


cell_filter_mod = cell_filter.copy() # this is a bool arr, so can add up

cell_filter_mod[:, 31] = cell_filter_mod[:, 31] + cell_filter_mod[:, 30] # pretend cells who prefer stim #31 also prefer stim #32, take union
cell_filter_mod[:, 30] = 0 # remove stim #30
# cell_filter_mod[:, 37] = cell_filter_mod[:, 37] + cell_filter_mod[:, 38]
# cell_filter_mod[:, 38] = 0
# cell_filter_mod[:, 37] = cell_filter_mod[:, 37] + cell_filter_mod[:, 39]
# cell_filter_mod[:, 39] = 0
cell_filter_mod[:, 38] = cell_filter_mod[:, 38] + cell_filter_mod[:, 39]
cell_filter_mod[:, 39] = 0

cell_filter_mod[:, 41] = cell_filter_mod[:, 41] + cell_filter_mod[:, 40]
cell_filter_mod[:, 40] = 0
# cell_filter_mod[:, 47] = cell_filter_mod[:, 47] + cell_filter_mod[:, 48]
# cell_filter_mod[:, 48] = 0
# cell_filter_mod[:, 47] = cell_filter_mod[:, 47] + cell_filter_mod[:, 49]
# cell_filter_mod[:, 49] = 0
cell_filter_mod[:, 48] = cell_filter_mod[:, 48] + cell_filter_mod[:, 49]
cell_filter_mod[:, 49] = 0

cell_filter_mod[cell_filter_mod>1] = 1 # binarize cell_filter_mod

print(np.sum(cell_filter[:, 30:40], axis=0))
print(np.sum(cell_filter[:, 40:50], axis=0))
print(np.sum(cell_filter_mod[:, 30:40], axis=0))
print(np.sum(cell_filter_mod[:, 40:50], axis=0))


# In[258]:


resp_ad_by_stim, resp_tg_by_stim = resp_by_stim(resp_ad, resp_tg, cell_filter_mod, stim_id_mod)
adp_stim_boot = calc_adp_pop_boot_by_stim(resp_ad_by_stim, nboot=1000)
print(adp_stim_boot.shape, resp_ad_by_stim[0].shape, resp_ad_by_stim.__len__())

spatial_freq_cpd = np.round(np.geomspace(0.03, 0.9, num=10), 2)
spatial_freq_cpd = spatial_freq_cpd[1:-1] # remove lowest & 2 highest SF
noise_nbyn = [2,3,5,17, 6,10,34, 15,51,85] # from do_image_preprocess.py -> 2x2 grid to 85x85 grid
noise_nbyn = noise_nbyn[1:-1]
cell_filter_mod = np.delete(cell_filter_mod, [30, 39, 40, 49], axis=1) # remove merged stim

adp_stim_boot_pop = np.nanmean(adp_stim_boot, axis=2) # [nboot x nstim], avg over cells
adp_stim_pop_mean = np.nanmean(adp_stim_boot_pop, axis=0) # [nstim], avg over boot
adp_stim_pop_err = np.nanstd(adp_stim_boot_pop, axis=0)
adp_stim_pop_CI95 = np.nanpercentile(adp_stim_boot_pop, [2.5, 97.5], axis=0) # [2 x nstim], 95% CI


# In[259]:


plot_adp_by_SF(cell_filter_mod, adp_stim_pop_mean, adp_stim_pop_CI95, adp_stim_pop_err)


# In[260]:


def plot_adp_by_SF_grat(cell_filter, adp_stim_pop_mean, adp_stim_pop_CI95, adp_stim_pop_err):
    x = np.arange(len(adp_stim_pop_mean[30:38]))
    y1 = adp_stim_pop_mean[30:38]
    y2 = np.sum(cell_filter[:, 30:38], axis=0)

    fig, ax1 = plt.subplots(figsize=(6,5))
    # ax2 = ax1.twinx()

    ax1.plot(x, y1, 'b')
    ax1.fill_between(x, y1-adp_stim_pop_err[30:38], y1+adp_stim_pop_err[30:38], alpha=0.3, color='b');
    ax1.fill_between(x, adp_stim_pop_CI95[0, 30:38], adp_stim_pop_CI95[1, 30:38], alpha=0.3, color='orange');
    ax1.axhline(0, color='k', linestyle='-', alpha=0.3)

    # ax2.plot(x, y2, 'b-', alpha=0.3)
    # ax2.axvline(len(x)//2 - 0.5, color='k', linestyle='-.', alpha=0.3)

    fontsize = 16
    ax1.set_xlabel('spatial frequency (cpd)', fontsize=fontsize)
    ax1.set_xticks(x, list(spatial_freq_cpd), rotation=45);
    ax1.set_ylabel('adaptation index', fontsize=fontsize)
    # ax2.set_ylabel(f'filtered ncell: {cell_filter_name}', color='b')

    ax1.grid(False)
    # ax2.grid(False) # grid off

    # ax1 right and top spines off
    ax1.spines['right'].set_visible(False)
    ax1.spines['top'].set_visible(False)

    plt.tight_layout()
    plt.show()
    fig.savefig(f'adp_by_SF_grat_{cell_filter_name}.pdf')

plot_adp_by_SF_grat(cell_filter_mod, adp_stim_pop_mean, adp_stim_pop_CI95, adp_stim_pop_err)


# ##### boxplot for poster SFN

# In[192]:


adp_stim_boot_grat = np.mean(adp_stim_boot[:, 30:38, :], axis=0).T
cell_filter_mod_grat = cell_filter_mod[:, 30:38]
adp_stim_boot_grat.shape, cell_filter_mod_grat.shape


# In[205]:


df_adp_SF = pd.DataFrame()

# construct SF column
ncell = adp_stim_boot_grat.shape[0]
SF = list(spatial_freq_cpd) * ncell
SF = np.sort(np.array(SF))
assert len(SF) == ncell * len(spatial_freq_cpd)

# for each SF, loop thru all cells. so adp column is flattened adp_stim_boot_grat
adp_SF = adp_stim_boot_grat.flatten('F') # flatten by column (SF)
assert len(adp_SF) == len(SF)

# cell filter column
cell_filter_SF = cell_filter_mod_grat.flatten('F')
assert len(cell_filter_SF) == len(SF)

df_adp_SF['SF'] = SF
df_adp_SF['adp'] = adp_SF
df_adp_SF['cell_filter'] = cell_filter_SF
df_adp_SF

test = df_adp_SF.dropna(inplace=False)
test['cell_filter'].unique()


# In[210]:


f, axes = plt.subplots(1, 2, figsize=(12, 6))

# sns.boxplot(
#     data=df, x="SF", y="adp",
#     notch=True, showcaps=False,
#     bootstrap=10000,
#     medianprops={"color": "gold"},
#     boxprops={"facecolor": '#78acbc'},
#     ax=axes[0],
# )
# axes[0].set_xticklabels(spatial_freq_cpd, rotation=45)
# axes[0].set_xlabel('Spatial frequency (cpd)')
# axes[0].set_ylabel('Adaptation index')
# axes[0].set_ylim(-1., 2.5)

sns.boxplot(
    data=df_adp_SF, x="SF", y="adp",
    notch=True, showcaps=False,
    bootstrap=10000,
    medianprops={"color": "gold"},
    boxprops={"facecolor": '#78acbc'},
    ax=axes[1],
)
axes[1].set_xticklabels(spatial_freq_cpd, rotation=45)
axes[1].set_yticklabels([])
axes[1].set_xlabel('Preferred spatial frequency (cpd)')
axes[1].set_ylabel('')
# axes[1].set_ylim(-1., 2.5)

plt.tight_layout()
plt.show()

# f.savefig('adaptation_vs_spatial_freq_less_color.pdf', bbox_inches='tight')


# In[139]:


# determine ncell or ntrial bottleneck
    # ncell is important, bc pop adp stability is much better than single cell adp
    # but ntrial is bottleneck, bc pop adp stability increases as ntrial_rep increase
    # should still record yuansi suggested grat_SF5 data, w only 5 gratings, each 300 reps


# ### boot adp for single cell
# #### sort cell by adp

# In[117]:


adp_stim_boot.shape # nboot x nstim x ncell
adp_stim_cell = np.nanmean(adp_stim_boot, axis=0) # [nstim x ncell], avg over boot
# plt.matshow(adp_stim_cell);
# plt.colorbar();

# sort cells by overall adp across stims
adp_stim_cell_mean = np.nanmedian(adp_stim_cell, axis=0) # [ncell], avg over stims
adp_stim_cell_order = np.argsort(adp_stim_cell_mean)[::-1] # sort in descending order

plt.figure(figsize=(10,10))
plt.subplot(2,1,1)
plt.plot(np.sort(adp_stim_cell_mean)[::-1]);
plt.xlabel('cell #')
plt.ylabel('median adp across stims')

plt.subplot(2,1,2)
im = plt.imshow(adp_stim_cell[:, adp_stim_cell_order])
plt.grid(False)
plt.colorbar(fraction=0.015, pad=0.04)

plt.xlabel('cell sorted by median adp across stims')
plt.ylabel('stim')
plt.xticks([]) # no axis ticks
plt.yticks([])
plt.tight_layout()


# # resp

# In[118]:


resp_ad_by_stim[0].shape, len(resp_ad_by_stim) # ncell x ntrial_rep x nstim
resp_ad_stim_cell = np.array([i.sum(axis=1) for i in resp_ad_by_stim]) # [ncell x nstim], sum over trial reps
resp_ad_stim_cell.shape # nstim x ncell

# plt.matshow(resp_ad_stim_cell)
# plt.grid(False)
# plt.xlabel('cell')
# plt.ylabel('stim')
# plt.colorbar();


# ### pref stim 
# take only cells preferring stim  
# assume cells prefer 1 natural img, 1 grat, 1 noise (1 preferred stim among each stim type)  
# using resp of img driven cells only

# In[77]:


# print(resp_ad_stim_cell.shape) # nstim x ncell

R1_by_type = {}
R1_by_type['nat'] = resp_ad_stim_cell[:30]
R1_by_type['grat'] = resp_ad_stim_cell[30:40]
R1_by_type['noise'] = resp_ad_stim_cell[40:]
# print(R1_by_type['nat'].shape, R1_by_type['grat'].shape, R1_by_type['noise'].shape)

pref_stim_by_type = {}
stim_offset = 0
for itype, resp_matrix in R1_by_type.items():
    
    # pref_stim = np.argmax(resp_matrix, axis=0) + stim_offset # for each cell, get peak resp stim. can't nanargmax due to All-NaN slice
    pref_stim = np.ones((resp_matrix.shape[1], 1)) * np.pi
    for icell in np.arange(resp_matrix.shape[1]):
        resp_cell = resp_matrix[:, icell]
        try:
            pref_stim[icell] = np.nanargmax(resp_cell) + stim_offset
        except:
            pref_stim[icell] = np.nan # if cell not responsive to any stim in this stim type, set pref stim to nan
    pref_stim = [elem[0] for elem in pref_stim]
    pref_stim = np.array(pref_stim)

    # print(pref_stim.shape) # ncell
    pref_stim_by_type[itype] = pref_stim
    stim_offset = stim_offset + resp_matrix.shape[0]


# In[78]:


# how many cells prefer each stim?
for itype, pref_stim in pref_stim_by_type.items():
    print(itype)
    print(np.unique(pref_stim, return_counts=True)) # count elements in pref_stim

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15,5), )

for itype, pref_stim in pref_stim_by_type.items():
    ax1.plot(np.unique(pref_stim, return_counts=True)[1][:-1], label=itype, linestyle='--')
ax1.legend(frameon=False)

ax2.bar(np.arange(3), [55, 82/3, 86/3], alpha=0.5)
# ax2.legend(frameon=False)


# In[79]:


# pref_filter: which cell prefer which stim
nstim = resp_ad_stim_cell.shape[0]
ncell = resp_ad_stim_cell.shape[1]
pref_filter = np.zeros((nstim + 1, ncell)) # padding 1 row at the bottom for nan

for itype, pref_stim in pref_stim_by_type.items():
    pref_stim[np.isnan(pref_stim)] = nstim # set nan to last row
    pref_stim_int = [int(elem) for elem in pref_stim]
    pref_filter[pref_stim_int, np.arange(ncell)] = 1

pref_filter = pref_filter[:-1, :] # remove padding row. same format as img_driven: bool mat of nstim x ncell
# plt.matshow(pref_filter);
# np.save(os.path.join(data_dir, 'pref_filter.npy'), pref_filter)

# how many cells prefer each stim?
ncell_pref_stim = pref_filter.sum(axis=1) # nstim
plt.bar(np.arange(30), ncell_pref_stim[:30]);
plt.bar(np.arange(30, 40), ncell_pref_stim[30:40]);
plt.bar(np.arange(40, 50), ncell_pref_stim[40:]);


# ## sort cells by overall resp across stims

# In[80]:


# sort cells by overall resp across stims

resp_ad_stim_cell_agg = np.nanmedian(resp_ad_stim_cell, axis=0) # [ncell], avg over stims
resp_ad_stim_cell_order = np.argsort(resp_ad_stim_cell_agg)[::-1] # sort in descending order

plt.figure(figsize=(10,10))
plt.subplot(2,1,1)
plt.plot(np.sort(resp_ad_stim_cell_agg)[::-1]);
plt.xlabel('cell #')
plt.ylabel('median resp across stims')

plt.subplot(2,1,2)
plt.imshow(resp_ad_stim_cell[:, resp_ad_stim_cell_order])
plt.grid(False)
plt.colorbar(fraction=0.015, pad=0.04)

plt.xlabel('cell sorted by median resp across stims')
plt.ylabel('stim')
plt.xticks([]) # no axis ticks
plt.yticks([])
plt.tight_layout()


# ## sort stims by overall resp across cells

# In[81]:


# sort stims by overall resp across stims

resp_ad_stim_cell_agg = np.nanmedian(resp_ad_stim_cell, axis=1) # [nstim], avg over cells
resp_ad_stim_cell_order = np.argsort(resp_ad_stim_cell_agg)[::-1] # sort in descending order
resp_ad_stim_cell_order

plt.plot(np.sort(resp_ad_stim_cell_agg)[::-1]);
plt.xlabel('stim sorted by pop resp')
plt.ylabel('median resp across cells')

plt.matshow(resp_ad_stim_cell[resp_ad_stim_cell_order, :])
plt.grid(False)
plt.colorbar()
plt.xlabel('cell')
plt.ylabel('stim sorted by pop resp')
plt.xticks([]) # no axis ticks
plt.yticks([])
# plt.tight_layout()


# ## sort cells by tuning curve: vis filter only
# there is no good tuning curve for mix50 stim

# In[82]:


resp_ad_by_stim, resp_tg_by_stim = resp_by_stim(resp_ad, resp_tg, cell_filter, stim_id)
resp_ad_by_stim[0].shape, len(resp_ad_by_stim) # ncell x ntrial_rep x nstim

resp_ad_stim_cell = np.array([i.sum(axis=1) for i in resp_ad_by_stim]) # [ncell x nstim], sum over trial reps
resp_ad_stim_cell.shape # nstim x ncell

resp_stim = np.nanmedian(resp_ad_stim_cell, axis=1) # avg over cells
plt.plot(resp_stim);

ncell = resp_ad_stim_cell.shape[1]
nstim = resp_ad_stim_cell.shape[0]
resp_ad_sorted = resp_ad_stim_cell[resp_ad_stim_cell_order, :]
# TODO: fix this when img_driven_only == True


# In[102]:


# for icell in np.arange(ncell):
#     plt.plot(resp_ad_sorted[:, icell], label=icell, alpha=0.9)
# # overplotting


# In[105]:


# # single cell tuning curve
# nrow = 10
# ncol = 10
# nplot = nrow * ncol
# nfig = int(np.ceil(ncell / nplot))
# data = resp_ad_stim_cell # resp_ad_sorted

# for ifig in np.arange(nfig):
#     plt.figure(figsize=(20,10))
#     for iplot in np.arange(nplot):
#         icell = ifig * nplot + iplot
#         if icell < ncell:
#             plt.subplot(nrow, ncol, iplot+1)
#             x = np.arange(len(data[:, icell]))
#             ymin = np.nanmin(data[:, icell])
#             ymax = np.nanmax(data[:, icell])
#             plt.plot(x, data[:, icell], label=icell, alpha=0.9) # plt.plot for img driven filter. scatter for pref filter
#             plt.xticks([])
#             plt.yticks([])
#             plt.xlim([0, nstim])
#             try:
#                 plt.ylim([ymin*0.9, ymax*1.1])
#             except:
#                 pass
#             plt.grid(False)
#             # plt.xlabel('stim')
#             # plt.ylabel('resp')
#             plt.title('cell %d' % icell)
#     plt.tight_layout()


# ### corr btw grat vs noise SF tuning curve

# In[83]:


tuning_corr = []
for icell in np.arange(ncell):
    grat_tuning = resp_ad_stim_cell[30:40, icell]
    noise_tuning = resp_ad_stim_cell[40:50, icell]
    tuning_corr_cell = np.corrcoef(grat_tuning, noise_tuning)[0,1]
    plt.plot(grat_tuning, noise_tuning, '.', label=icell, alpha=0.4)
    tuning_corr.append(tuning_corr_cell)

plt.xlabel('grat tuning')
plt.ylabel('noise tuning')
plt.xscale('symlog')
plt.yscale('symlog')
plt.axis('equal');


# In[85]:


plt.hist(tuning_corr, bins=50);
plt.xlim(-1, 1);


# In[88]:


# sort cells by tuning_corr
tuning_corr = np.array(tuning_corr)
tuning_corr_order = np.argsort(tuning_corr)[::-1] # sort in descending order

plt.figure(figsize=(10,5))
plt.subplot(1,2,1)
plt.plot(resp_ad_stim_cell[30:40, tuning_corr_order[0]], label='grat tuning')
plt.plot(resp_ad_stim_cell[40:, tuning_corr_order[0]], label='noise tuning')
plt.text(0.7, 0.6, 'corr = %.2f' % tuning_corr[tuning_corr_order[0]], transform=plt.gca().transAxes, fontsize=12)
plt.legend(frameon=False)

plt.subplot(1,2,2)
plt.plot(resp_ad_stim_cell[30:40, tuning_corr_order[-1]], label='grat tuning')
plt.plot(resp_ad_stim_cell[40:, tuning_corr_order[-1]], label='noise tuning')
plt.text(0.7, 0.6, 'corr = %.2f' % tuning_corr[tuning_corr_order[-1]], transform=plt.gca().transAxes, fontsize=12)
plt.legend(frameon=False)

plt.tight_layout()


# In[137]:


# # single cell tuning curve corr btw grat and noise
# nrow = 10
# ncol = 10
# nplot = nrow * ncol
# nfig = int(np.ceil(ncell / nplot))

# tuning_corr_sort = tuning_corr[tuning_corr_order]
# resp_sortby_SF_corr = resp_ad_stim_cell[:, tuning_corr_order]
# data = resp_sortby_SF_corr

# for ifig in np.arange(nfig):
#     plt.figure(figsize=(20,10))
#     for iplot in np.arange(nplot):
#         icell = ifig * nplot + iplot
#         if icell < ncell:
#             plt.subplot(nrow, ncol, iplot+1)
#             x = np.arange(len(data[:, icell]))
#             plt.plot(resp_ad_stim_cell[30:40, icell], label='grat tuning')
#             plt.plot(resp_ad_stim_cell[40:, icell], label='noise tuning')
#             plt.xticks([])
#             plt.grid(False)
#             plt.title('cell %d' % icell + ' corr = %.2f' % tuning_corr_sort[icell])
#     plt.tight_layout()

