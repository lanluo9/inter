#!/usr/bin/env python
# coding: utf-8

# # prep

# In[1]:


import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
import plotly.express as px

import os
import pandas as pd
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


def load_data(dir_path):
    file_path = dir_path.replace('\\', '/')
    data = sio.loadmat(file_path + '/trace_trial_stim.mat')
    
    stim_seq = data['stim_seq']
    stim_id = [i[0] for i in stim_seq]
    trace_by_trial = data['trace_by_trial']

    ncell = trace_by_trial.shape[0]
    nstim = len(np.unique(stim_id))
    ntrial = trace_by_trial.shape[1]
    nframe = trace_by_trial.shape[2]
    print(ncell, nstim, ntrial, nframe)

    return stim_id, trace_by_trial, # ncell, nstim, ntrial, nframe


def calc_trace_stim(trace_by_trial, stim_id):
    trace_avg_cell = np.mean(np.mean(trace_by_trial, axis=0), axis=0)
    trace_stim_avg = []
    # trace_stim_std = []
    # trace_stim_sem = []

    for i in np.unique(stim_id):
        trace_istim_avg = np.mean(trace_by_trial[:, np.where(stim_id == i)[0]], axis=1) # ncell x nframe
        trace_istim_avg = np.mean(trace_istim_avg, axis=0) # nframe
        # trace_istim_std = np.std(trace_by_trial[:, np.where(stim_id == i)[0]], axis=1)
        # trace_istim_sem = trace_istim_std / np.sqrt(len(np.where(stim_id == i)[0]))

        trace_stim_avg.append(trace_istim_avg)
        # trace_stim_std.append(trace_istim_std)
        # trace_stim_sem.append(trace_istim_sem)

    print(len(trace_stim_avg), trace_stim_avg[0].shape)
    return trace_avg_cell, trace_stim_avg

def load_dfof(iset):
    dir_sub = str(meta.area[iset]) + '_i' + str(meta.mouse[iset]) + '_' + str(meta.date[iset])
    dfof = sio.loadmat(os.path.join(dir_name, dir_sub, 'dfof' + '.mat'))
    dfof_ad = dfof['dfof_ad']
    dfof_tg = dfof['dfof_tg']
    ncell = dfof['dfof_ad'].shape[0]
    return dfof_ad, dfof_tg, ncell

def threshold_dfof(dfof_ad, dfof_tg, thres_perc_low=1, thres_perc_high=95):
    thres_low = np.percentile(dfof_ad.flatten(), thres_perc_low)
    thres_high = np.percentile(dfof_ad.flatten(), thres_perc_high)
    dfof_ad_thres = dfof_ad.copy()
    dfof_ad_thres[dfof_ad_thres < thres_low] = np.nan # turn off left thres bc right tail is heavy
    dfof_ad_thres[dfof_ad_thres > thres_high] = np.nan
    dfof_tg_thres = dfof_tg.copy()
    dfof_tg_thres[dfof_tg_thres < thres_low] = np.nan # threshold resp_tg by resp_tg
    dfof_tg_thres[dfof_tg_thres > thres_high] = np.nan
    return dfof_ad_thres, dfof_tg_thres

def calc_adp(dfof_ad, dfof_tg, thres=np.nan):
    adp = (dfof_tg - dfof_ad) / (dfof_tg + dfof_ad + 1e-7) # changed definition of adp to (tg - ad) / (tg + ad)
    if thres is not np.nan:
        adp[adp > thres] = np.nan
        adp[adp < -thres] = np.nan
    return adp


# In[3]:


def plot_dfof(dfof_ad_cellpose, dfof_tg_cellpose, dfof_ad_manual, dfof_tg_manual):
    dfof_ad_cellpose_median = np.nanmedian(dfof_ad_cellpose.flatten())
    dfof_ad_manual_median = np.nanmedian(dfof_ad_manual.flatten())
    dfof_tg_cellpose_median = np.nanmedian(dfof_tg_cellpose.flatten())
    dfof_tg_manual_median = np.nanmedian(dfof_tg_manual.flatten())

    # xlow = -100
    # xhigh = 200

    plt.figure(figsize=(12, 10))
    plt.subplot(2,2,1)
    plt.hist(dfof_ad_cellpose.flatten(), bins=100, alpha=0.5, density=True, color='b');
    plt.axvline(x=dfof_ad_cellpose_median, color='b', linestyle='--')
    plt.hist(dfof_ad_manual.flatten(), bins=100, alpha=0.5, density=True, color='cyan');
    plt.axvline(x=dfof_ad_manual_median, color='cyan', linestyle='--')
    # plt.xlim([xlow, xhigh])
    plt.legend(['cellpose ad', 'manual ad'])

    plt.subplot(2,2,2)
    plt.hist(dfof_tg_cellpose.flatten(), bins=100, alpha=0.5, density=True, color='orange');
    plt.axvline(x=dfof_tg_cellpose_median, color='orange', linestyle='--')
    plt.hist(dfof_tg_manual.flatten(), bins=100, alpha=0.5, density=True, color='yellow');
    plt.axvline(x=dfof_tg_manual_median, color='yellow', linestyle='--')
    # plt.xlim([xlow, xhigh])
    plt.legend(['cellpose tg', 'manual tg'])

    plt.subplot(2,2,3)
    plt.hist(dfof_ad_cellpose.flatten(), bins=100, alpha=0.5, density=True, color='b');
    plt.axvline(x=dfof_ad_cellpose_median, color='b', linestyle='--')
    plt.hist(dfof_tg_cellpose.flatten(), bins=100, alpha=0.5, density=True, color='orange');
    plt.axvline(x=dfof_tg_cellpose_median, color='orange', linestyle='--')
    # plt.xlim([xlow, xhigh])
    plt.legend(['cellpose ad', 'cellpose tg'])

    plt.subplot(2,2,4)
    plt.hist(dfof_ad_manual.flatten(), bins=100, alpha=0.5, density=True, color='cyan');
    plt.axvline(x=dfof_ad_manual_median, color='cyan', linestyle='--')
    plt.hist(dfof_tg_manual.flatten(), bins=100, alpha=0.5, density=True, color='yellow');
    plt.axvline(x=dfof_tg_manual_median, color='yellow', linestyle='--')
    # plt.xlim([xlow, xhigh])
    plt.legend(['manual ad', 'manual tg']);


def plot_adp_compare(adp_cellpose, adp_manual):
    plt.figure(figsize=(15, 6))

    plt.subplot(1,2,1)
    plt.plot(np.nanmean(adp_cellpose, axis=0), color='b', alpha=0.5)
    plt.plot(np.nanmean(adp_manual, axis=0), color='orange')
    plt.axhline(y=0, color='gray', linestyle='--')
    plt.axhline(y=np.nanmean(adp_cellpose.flatten()), color='b', linestyle='--', alpha=0.5)
    plt.axhline(y=np.nanmean(adp_manual.flatten()), color='orange', linestyle='--')
    plt.legend(['cellpose adp mean', 'manual adp mean'], loc='lower left')

    plt.subplot(1,2,2)
    plt.plot(np.nanmedian(adp_cellpose, axis=0), color='b', alpha=0.5)
    plt.plot(np.nanmedian(adp_manual, axis=0), color='orange')
    plt.axhline(y=0, color='gray', linestyle='--')
    plt.axhline(y=np.nanmedian(adp_cellpose.flatten()), color='b', linestyle='--', alpha=0.5)
    plt.axhline(y=np.nanmedian(adp_manual.flatten()), color='orange', linestyle='--')
    # plt.ylim([-0.5,0.25])
    plt.legend(['cellpose adp median', 'manual adp median'], loc='lower left');


# # compare trace

# In[4]:


# i1369 V1 high res, cellpose
dir_path = r'Z:\All_Staff\home\lan\Data\2P_images\mat_inter\V1_i1369_220310_cellpose'
stim_id_cellpose, trace_by_trial_cellpose = load_data(dir_path)
trace_avg_cell_cellpose, trace_stim_avg_cellpose = calc_trace_stim(trace_by_trial_cellpose, stim_id_cellpose)

# i1369 V1 high res, manual
dir_path = r'Z:\All_Staff\home\lan\Data\2P_images\mat_inter\V1_i1369_220310'
stim_id_manual, trace_by_trial_manual = load_data(dir_path)
trace_avg_cell_manual, trace_stim_avg_manual = calc_trace_stim(trace_by_trial_manual, stim_id_manual)

# normalize trace to same max min
trace_avg_cell_cellpose_norm = (trace_avg_cell_cellpose - np.min(trace_avg_cell_cellpose)) / (np.max(trace_avg_cell_cellpose) - np.min(trace_avg_cell_cellpose))
trace_avg_cell_manual_norm = (trace_avg_cell_manual - np.min(trace_avg_cell_manual)) / (np.max(trace_avg_cell_manual) - np.min(trace_avg_cell_manual))


# In[19]:


plt.figure(figsize=(12, 5))

# plt.subplot(2,1,1) # no normalization
# plt.plot(trace_avg_cell_cellpose, alpha=0.8, linewidth=3) # cellpose higher by 50
# plt.plot(trace_avg_cell_manual, alpha=0.8, linewidth=3)

# plt.axvline(x=np.argmax(trace_avg_cell_cellpose[:25]), color='gray', linestyle='--')
# plt.axvline(x=np.argmin(trace_avg_cell_cellpose[25:30]) + 25, color='gray', linestyle='--')
# plt.axvline(x=np.argmax(trace_avg_cell_cellpose), color='gray', linestyle='--')
# plt.legend(['cellpose', 'manual'])

# plt.subplot(2,1,2) # normalized via maxmin
plt.plot(trace_avg_cell_cellpose_norm, alpha=0.8, linewidth=3, linestyle='--')
plt.plot(trace_avg_cell_manual_norm, alpha=0.5, linewidth=3)
plt.axvline(11, color='gray', linestyle='--')
plt.axvline(13, color='gray', linestyle='--')
plt.legend(['cellpose norm', 'manual norm']);


# # compare resp

# In[8]:


metadata = {'mouse': [1369, 1369], 'date': ['220310', '220310_cellpose'], 'area': ['V1', 'V1']}
meta = pd.DataFrame(data=metadata)
print(meta)

nset = len(meta.index); ncell = []; nori = 30; nframe_trial = 156
dir_name = r'Z:\All_Staff\home\lan\Data\2P_images\mat_inter'.replace('\\', '/')


# ## resp dist: not gauss or log gauss

# In[27]:


resp_cells = trace_by_trial_cellpose[:,:,11:14].mean(axis=2)
resp_pop = resp_cells.mean(axis=0)
# resp_pop.shape
plt.hist(resp_pop, bins=100, alpha=0.5);


# In[32]:


from scipy import stats
print(stats.kstest(resp_pop, 'norm')) # p value almost 0, so not gaussian

resp_pop_pos = resp_pop - np.min(resp_pop) + 1e-7
resp_pop_log = np.log(resp_pop_pos)
stats.kstest(resp_pop_log, 'norm')


# In[33]:


ncell = resp_cells.shape[0]
for i in range(ncell):
    resp_cell = resp_cells[i]
    # resp_cell_pos = resp_cell - np.min(resp_cell) + 1e-7
    # resp_cell_log = np.log(resp_cell_pos)

    st_cell, p_cell = stats.kstest(resp_cell, 'norm')
    if p_cell > 0.001:
        print(st_cell, p_cell)


# ## extreme resp rarely co-occur
# check if extreme value of dfof_ad and dfof_tg co-occur -> not really  
# conclusion: should threshold both dfof_ad and dfof_tg at top / bottom 5% (substitute with np.nan)

# In[222]:


dfof_ad_cellpose, dfof_tg_cellpose, ncell_cellpose = load_dfof(1)
dfof_ad_manual, dfof_tg_manual, ncell_manual = load_dfof(0)
large_ad_cellpose_id = dfof_ad_cellpose.flatten() > np.percentile(dfof_ad_cellpose.flatten(), 99)
large_ad_manual_id = dfof_ad_manual.flatten() > np.percentile(dfof_ad_manual.flatten(), 99)

large_ad_cellpose = dfof_ad_cellpose.flatten()[large_ad_cellpose_id]
large_ad_manual = dfof_ad_manual.flatten()[large_ad_manual_id]
large_ad_tg_cellpose = dfof_tg_cellpose.flatten()[large_ad_cellpose_id] # tg corresponding to large_ad
large_ad_tg_manual = dfof_tg_manual.flatten()[large_ad_manual_id]

plt.hist(large_ad_cellpose.flatten(), bins=30, alpha=0.5, density=True, color='b');
plt.hist(large_ad_manual.flatten(), bins=30, alpha=0.5, density=True, color='cyan');
plt.hist(large_ad_tg_cellpose.flatten(), bins=30, alpha=0.5, density=True, color='orange');
plt.hist(large_ad_tg_manual.flatten(), bins=30, alpha=0.5, density=True, color='yellow');


# ## resp dist

# In[282]:


# check distribution of response to adapter and target, original value w/o thresholding

dfof_ad_cellpose, dfof_tg_cellpose, ncell_cellpose = load_dfof(1)
dfof_ad_manual, dfof_tg_manual, ncell_manual = load_dfof(0)
plot_dfof(dfof_ad_cellpose, dfof_tg_cellpose, dfof_ad_manual, dfof_tg_manual)


# ## find resp thres

# In[287]:


thres_perc_low = 1
thres_perc_high = 90

dfof_ad_manual, dfof_tg_manual, ncell_manual = load_dfof(0)
dfof_ad_manual, dfof_tg_manual = threshold_dfof(dfof_ad_manual, dfof_tg_manual, 
                                                thres_perc_low=thres_perc_low, thres_perc_high=thres_perc_high)

dfof_ad_cellpose, dfof_tg_cellpose, ncell_cellpose = load_dfof(1)
dfof_ad_cellpose, dfof_tg_cellpose = threshold_dfof(dfof_ad_cellpose, dfof_tg_cellpose, 
                                                    thres_perc_low=thres_perc_low, thres_perc_high=thres_perc_high)

plot_dfof(dfof_ad_cellpose, dfof_tg_cellpose, dfof_ad_manual, dfof_tg_manual)


# # compare adp

# In[290]:


dfof_ad_cellpose, dfof_tg_cellpose, ncell_cellpose = load_dfof(1)
dfof_ad_cellpose, dfof_tg_cellpose = threshold_dfof(dfof_ad_cellpose, dfof_tg_cellpose, thres_perc_low=thres_perc_low, thres_perc_high=thres_perc_high)
adp_cellpose = calc_adp(dfof_ad_cellpose, dfof_tg_cellpose, thres=5) # threshold adp to rid of extreme values

dfof_ad_manual, dfof_tg_manual, ncell_manual = load_dfof(0)
dfof_ad_manual, dfof_tg_manual = threshold_dfof(dfof_ad_manual, dfof_tg_manual, thres_perc_low=thres_perc_low, thres_perc_high=thres_perc_high)
adp_manual = calc_adp(dfof_ad_manual, dfof_tg_manual, thres=5)

print(adp_cellpose.shape, adp_manual.shape)
np.nanmedian(adp_manual.flatten()), np.nanmedian(adp_cellpose.flatten()), np.nanmean(adp_manual.flatten()), np.nanmean(adp_cellpose.flatten())


# In[291]:


plt.figure(figsize=(12, 5))
plt.hist(adp_cellpose.flatten(), bins=100, alpha=0.5, density=True, color='b');
plt.axvline(x=np.nanmedian(adp_cellpose.flatten()), color='b', linestyle='--')
print(f'cellpose adp median: {np.nanmedian(adp_cellpose.flatten()):.2f}, cellpose mean: {np.nanmean(adp_cellpose.flatten()):.2f}')

plt.hist(adp_manual.flatten(), bins=100, alpha=0.5, density=True, color='orange');
plt.axvline(x=np.nanmedian(adp_manual.flatten()), color='orange', linestyle='--')
print(f'manual adp median: {np.nanmedian(adp_manual.flatten()):.2f}, manual mean: {np.nanmean(adp_manual.flatten()):.2f}')

plt.legend(['cellpose adp', 'manual adp']);
# plt.yscale('log') # use if density=False
# plt.xlim([-200, 200]); # use if adp not thresholded


# In[292]:


# df1 = pd.DataFrame(data=adp_cellpose.flatten(), columns=['adp'])
# df1['tag'] = 'cellpose'
# df2 = pd.DataFrame(data=adp_manual.flatten(), columns=['adp'])
# df2['tag'] = 'manual'

# df = pd.concat([df1, df2], axis=0)
# ax = sns.boxplot(x="adp", y="tag", data=df) # , whis=np.inf


# ## adp thres is necessary
# thresholding resp_ad and tg does not help
# ### 3d scatter: resp ad, resp tg, adp

# In[340]:


dfof_ad_cellpose, dfof_tg_cellpose, ncell_cellpose = load_dfof(1)
dfof_ad_cellpose, dfof_tg_cellpose = threshold_dfof(dfof_ad_cellpose, dfof_tg_cellpose, thres_perc_low=thres_perc_low, thres_perc_high=thres_perc_high)
adp_cellpose = calc_adp(dfof_ad_cellpose, dfof_tg_cellpose) # threshold adp to rid of extreme values

dfof_ad_manual, dfof_tg_manual, ncell_manual = load_dfof(0)
dfof_ad_manual, dfof_tg_manual = threshold_dfof(dfof_ad_manual, dfof_tg_manual, thres_perc_low=thres_perc_low, thres_perc_high=thres_perc_high)
adp_manual = calc_adp(dfof_ad_manual, dfof_tg_manual)

fig = plt.figure(figsize=(12, 5))
ax = fig.add_subplot(111, projection='3d')
ax.scatter(dfof_ad_cellpose.flatten(), dfof_tg_cellpose.flatten(), adp_cellpose.flatten(), c='b', marker='o', alpha=0.2, s=10)
ax.scatter(dfof_ad_manual.flatten(), dfof_tg_manual.flatten(), adp_manual.flatten(), c='orange', marker='o', alpha=0.8, s=10)
ax.set_xlabel('dfof ad')
ax.set_ylabel('dfof tg')
ax.set_zlabel('adp');

ax.view_init(elev=15, azim=30, ) # set angle of view


# In[316]:


fig = px.scatter_3d(x=dfof_ad_cellpose.flatten(), y=dfof_tg_cellpose.flatten(), z=adp_cellpose.flatten(),)
fig = px.scatter_3d(x=dfof_ad_manual.flatten(), y=dfof_tg_manual.flatten(), z=adp_manual.flatten(),)
fig.update_traces(marker_size = 2)
fig.show()


# ## compare adp by stim

# In[341]:


# without thres adp, plot adp compare across stims

dfof_ad_cellpose, dfof_tg_cellpose, ncell_cellpose = load_dfof(1)
dfof_ad_cellpose, dfof_tg_cellpose = threshold_dfof(dfof_ad_cellpose, dfof_tg_cellpose)
adp_cellpose = calc_adp(dfof_ad_cellpose, dfof_tg_cellpose) # no threshold, see if groupby stim works

dfof_ad_manual, dfof_tg_manual, ncell_manual = load_dfof(0)
dfof_ad_manual, dfof_tg_manual = threshold_dfof(dfof_ad_manual, dfof_tg_manual)
adp_manual = calc_adp(dfof_ad_manual, dfof_tg_manual)

print(adp_cellpose.shape, adp_manual.shape)
print(np.nanmedian(adp_manual.flatten()), np.nanmedian(adp_cellpose.flatten()), np.nanmean(adp_manual.flatten()), np.nanmean(adp_cellpose.flatten()))
plot_adp_compare(adp_cellpose, adp_manual)


# In[342]:


# without thres adp, plot adp compare across stims

dfof_ad_cellpose, dfof_tg_cellpose, ncell_cellpose = load_dfof(1)
dfof_ad_cellpose, dfof_tg_cellpose = threshold_dfof(dfof_ad_cellpose, dfof_tg_cellpose) # cellpose did not filter by vis-driven
adp_cellpose = calc_adp(dfof_ad_cellpose, dfof_tg_cellpose, thres=5) # adp thres is still necessary even if groupby stim

dfof_ad_manual, dfof_tg_manual, ncell_manual = load_dfof(0)
dfof_ad_manual, dfof_tg_manual = threshold_dfof(dfof_ad_manual, dfof_tg_manual)
adp_manual = calc_adp(dfof_ad_manual, dfof_tg_manual, thres=5)

print(adp_cellpose.shape, adp_manual.shape)
print(np.nanmedian(adp_manual.flatten()), np.nanmedian(adp_cellpose.flatten()), np.nanmean(adp_manual.flatten()), np.nanmean(adp_cellpose.flatten()))
plot_adp_compare(adp_cellpose, adp_manual)

