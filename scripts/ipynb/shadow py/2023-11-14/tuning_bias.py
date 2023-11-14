#!/usr/bin/env python
# coding: utf-8

# # prep

# In[33]:


## on nuke, we use env `anaconda3` due to `base` being outdated

import numpy as np
from numpy import dot
from numpy.linalg import norm

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

import scipy.io as sio
import scipy.stats as stats
from statsmodels.stats.proportion import proportion_confint

from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.svm import LinearSVC
from sklearn.model_selection import train_test_split
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

from tqdm import tqdm
from IPython.display import clear_output
import os
import pickle

get_ipython().run_line_magic('load_ext', 'autoreload')
get_ipython().run_line_magic('autoreload', '2')


# In[2]:


local_flag = False
if local_flag:
    repo_dir = r'D:\repo\inter_data\inter'.replace("\\", "/") # under env dimred
else:
    # repo_dir = r'C:\Users\ll357\Documents\inter'.replace("\\", "/")
    repo_dir = r'C:\Users\lan\Documents\repos\inter'.replace("\\", "/")
os.chdir(repo_dir)
from src import adp


# In[4]:


dir_inter = r'\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images\mat_inter/'.replace('\\', '/')
dir_file = os.path.join(dir_inter, 'adp_dataset_master.xlsx')
data_info = pd.read_excel(dir_file)
data_info.head()

meta = data_info[(data_info.paradigm == 'grating') 
                 & (data_info.area == 'LM')
                 & ((data_info.cellpose_seg == True) | (data_info.manual_seg == True))] # ensure segmentation
meta = meta.reset_index(drop=True)
nset = meta.shape[0]
print(meta.area.value_counts(), nset)

meta_LM = meta.copy()
meta_LM.tail()


# In[5]:


meta = data_info[(data_info.paradigm == 'grating') 
                 & (data_info.area == 'V1') 
                 & (data_info.gcamp == '6s') # avoid mixing in gcamp8f
                 & ((data_info.cellpose_seg == True) | (data_info.manual_seg == True))] # there are 2 data in V1 that haven't been segmented due to cellpose failing with low quality tiff
meta = meta.reset_index(drop=True)
nset = meta.shape[0]
print(meta.area.value_counts(), nset)

meta_V1 = meta.copy()
meta_V1


# In[6]:


meta = data_info[(data_info.paradigm == 'grating') 
                 & (data_info.area == 'LI') 
                 & (data_info.gcamp == '6s') # avoid mixing in gcamp8f
                 & (data_info.manual_seg != 'TODO') # 2 LI data still need manual segm
                 & (data_info.note.str.contains('bad') != True) # exclude bad data
                 ]
meta = meta.reset_index(drop=True)
nset = meta.shape[0]
print(meta.area.value_counts(), nset)

meta_LI = meta.copy()
meta_LI


# # batch write df_tidy (don't rerun unless needed)

# In[1059]:


meta = pd.concat([meta_V1, meta_LM, meta_LI], axis=0).reset_index(drop=True)
nset = len(meta)

for iset in tqdm(range(nset)):
    print(f'iset={iset}, nset={nset}')


    ## load data
    mouse = meta.loc[iset, 'mouse'].astype(str)
    imouse = 'i' + mouse
    date = meta.loc[iset, 'date'].astype(str)
    area = meta.loc[iset, 'area']
    sess = meta.loc[iset, 'num']
    print(mouse, date, area, sess)
    dir_identifier = f'{area}_{imouse}_{date}_{sess}'

    mat_inter = r'\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images\mat_inter/'.replace('\\', '/')
    for dirname in os.listdir(mat_inter):
        if dir_identifier in dirname:
            dir_data = mat_inter + dirname
            break

    stim_id, resp_ad, resp_tg = adp.load_resp_trial(os.path.join(dir_data), vis_filter=False) # already sliced by resp time window from matlab. only differentiating isi in resp_tg by indexing into trial. match isi_nframe
    ## NOTE: stim_id contains the final trial, which should be cut off due to incomplete trial, to match trace_by_trial


    ## construct dfof
    R1_dfof = resp_ad.flatten() # cell x trial, row major flatten
    R2_dfof = resp_tg.flatten()
    dfof = np.hstack((R1_dfof, R2_dfof)) # sort by cell_id, then trial_id, finally resp_id


    ## construct cell, trial, resp id
    ncell = resp_ad.shape[0]
    ntrial = resp_ad.shape[1]
    cell_id = np.repeat(np.arange(ncell), ntrial) # repeat cell id arr for ntrial times == sort by cell_id, then trial_id
    cell_id = np.hstack((cell_id, cell_id)) # stack two copies for R1 and R2. finally sort by resp_id

    trial_id = np.tile(np.arange(ntrial), ncell) # for each cell, go over all trials
    trial_id = np.hstack((trial_id, trial_id))

    len_df = ncell * ntrial * 2 # 2 for (R1, R2)
    resp_id = ['R1'] * int(len_df/2) + ['R2'] * int(len_df/2) # first half is flattened resp_ad, second half is flattened resp_tg


    ## construct stim info col: stim 2 orientation, stim 1 orien, isi, adapter contrast
    trial_stim_orien = [item[0] for item in stim_id['stim_ori']] # print(np.sort(np.unique(trial_stim_orien)))
    trial_stim_orien_dict = {} ## map stim2 orientation to int. convert from grat ori in deg to 0-based indexing (ensured)
    for i, item in enumerate(np.sort(np.unique(trial_stim_orien))):
        trial_stim_orien_dict[item] = i
    trial_stim_orien_int = [trial_stim_orien_dict[item] for item in trial_stim_orien]

    trial_isi_nframe = [item[0] for item in stim_id['isi_nframe']]
    trial_adapter_contrast = [item[0] for item in stim_id['adapter_contrast']]

    trial_stim_orien_int = trial_stim_orien_int[:ntrial] # if any stim info longer than ntrial, slice off the last one
    trial_isi_nframe = trial_isi_nframe[:ntrial]
    trial_adapter_contrast = trial_adapter_contrast[:ntrial]
    

    ## make stim info col: same as trial_id - tile then hstack
    stim_id_col = np.tile(trial_stim_orien_int, ncell)
    isi_col = np.tile(trial_isi_nframe, ncell)
    ad_con_col = np.tile(trial_adapter_contrast, ncell)

    stim_id_col = np.hstack((stim_id_col, stim_id_col)) # stim2 orientation
    isi_col = np.hstack((isi_col, isi_col))
    ad_con_col = np.hstack((ad_con_col, ad_con_col))

    df_tidy = pd.DataFrame({'dfof': dfof, 'cell_id': cell_id, 'trial_id': trial_id, 'resp_id': resp_id, 
                            'isi': isi_col, 'stim1_contrast': ad_con_col, 
                            'stim2_id': stim_id_col,})
    df_tidy['area'] = area
    df_tidy['stim1_id'] = 0 # adapter (stim1) is always 0 deg / vertical gratings
    df_tidy.loc[(df_tidy.stim1_contrast == 0) & (df_tidy.resp_id == 'R1'), 'dfof'] = np.nan # set R1 to nan if no adapter (stim1_contrast == 0)
    df_tidy = df_tidy.dropna().reset_index(drop=True) # drop nan rows

    df_tidy['isi'] = df_tidy['isi'].apply(lambda x: 250 if x <= 10 else 750) # convert isi_nframe to ms. 250 ms is about 8 frames, 750 ms is about 22-23 frames
    df_tidy.loc[df_tidy.stim1_contrast == 0, 'isi'] = 6000 # set ISI to 6000 ms (equal to ITI) if no adapter (stim1_contrast == 0)


    ## vis cell filter, well_fit filter, & img driven cell-stim filter
    df_tidy['filter_cell_vis'] = np.nan
    df_tidy['filter_cell_vis_pval'] = np.nan # allow continuous filtering on how significant the cell gets visually driven
    df_tidy['filter_cell_well_fit'] = np.nan
    df_tidy['filter_cell_stim'] = np.nan

    with open(os.path.join(dir_data, 'vis_driven_ttest_bonferroni.pickle'), 'rb') as f: # changed to strict bonferroni
    # with open(os.path.join(dir_data, 'vis_driven_ttest_bonferroni_strict.pickle'), 'rb') as f: # changed to strict bonferroni
        filter_file = pickle.load(f)
    filter_cell_stim = filter_file['img_driven']
    filter_cell_vis = filter_file['vis_driven']
    filter_cell_vis_pval = np.min(filter_file['p_ttest'], axis=1) # min pval across all stim
        
    well_fit = sio.loadmat(os.path.join(dir_data, 'fit_bootstrap_90perc.mat'))
    well_fit_cell = np.array([x[0] for x in well_fit['well_fit_cell']])

    for icell in np.arange(filter_cell_stim.shape[0]): # filter_cell_stim is ncell x nstim
        df_tidy.loc[df_tidy['cell_id']==icell, 'filter_cell_vis'] = filter_cell_vis[icell]
        df_tidy.loc[df_tidy['cell_id']==icell, 'filter_cell_vis_pval'] = filter_cell_vis_pval[icell]
        df_tidy.loc[df_tidy['cell_id']==icell, 'filter_cell_well_fit'] = well_fit_cell[icell]
        for istim in np.arange(filter_cell_stim.shape[1]):
            df_tidy.loc[(df_tidy['stim2_id']==istim) & (df_tidy['cell_id']==icell), 'filter_cell_stim'] = filter_cell_stim[icell, istim]
    # df_tidy.filter_cell_vis.value_counts(), df_tidy.filter_cell_stim.value_counts()


    # ## cell tuning in 3 possible ISI, 2 runs (each run uses half of the trials)
    # fit_tuning_half = sio.loadmat(os.path.join(dir_data, 'fit_tuning_half_trials.mat'))
    # ori_pref_run1 = fit_tuning_half['ori_pref_runs'][:, :, 0] # ncell x 3isi x 2run, [noad vs ad750 vs ad250]
    # ori_pref_run2 = fit_tuning_half['ori_pref_runs'][:, :, 1]

    # ori_pref_runs_sorted = well_fit['ori_pref_runs']
    # # ori_pref_runs_sorted = np.array([np.sort(x) for x in well_fit['ori_pref_runs']]) # ncell x nrun. sort each row of ori_pref_runs
    # # for icell in np.arange(10):
    # #     plt.plot(ori_pref_runs_sorted[icell, :])

    # percentile_threshold = 0.9
    # # if area == 'LI':
    # #     percentile_threshold = 0.7 # taken from well_fit_cell_criteria(_relax).m, only relax for LI
    # nrun = ori_pref_runs_sorted.shape[1]
    # rand_idx = np.random.randint(nrun*(1-percentile_threshold), nrun*percentile_threshold, size=2) # # bc well_fit cells are defined as: 90% of bootstrapped ori_pref within 22.5 deg, randomly select 2 indices between 100-900 in ori_pref_runs_sorted
    # ori_pref_noad1 = ori_pref_runs_sorted[:, rand_idx[0]] # ncell x 1
    # ori_pref_noad2 = ori_pref_runs_sorted[:, rand_idx[1]]


    # ## goodness of fit (R square) in isi 250 or 750
    # fit_tuning = sio.loadmat(os.path.join(dir_data, 'fit_tuning_isi3.mat')) # fit_tuning['fit_param'].shape # ncell x nparam x nstim
    # R_square = fit_tuning['fit_param'][:, -1, :] # ncell x nstim, final param is R_square of fit
    # R_square_750 = R_square[:, 1] # use R_sq to determine well_fit_ad
    # R_square_250 = R_square[:, 2]
    # well_fit_ad_250 = R_square_250 >= np.percentile(R_square_250, 50) # ncell x 1. only take top 50% of R_sq
    # well_fit_ad_750 = R_square_750 >= np.percentile(R_square_750, 50)


    # ## write cell property to df_tidy
    # ncell = df_tidy.cell_id.unique().shape[0]
    # for icell in np.arange(ncell):
    #     df_tidy.loc[(df_tidy.cell_id == icell), 'ori_pref_noad'] = ori_pref_run1[icell, 0] # ori_pref_run1 is ncell x 3isi, [noad vs ad750 vs ad250]
    #     df_tidy.loc[(df_tidy.cell_id == icell), 'ori_pref_ad_750'] = ori_pref_run1[icell, 1]
    #     df_tidy.loc[(df_tidy.cell_id == icell), 'ori_pref_ad_250'] = ori_pref_run1[icell, 2]

    #     df_tidy.loc[(df_tidy.cell_id == icell), 'ori_pref_noad1'] = ori_pref_run1[icell, 0] # for control plot. ori_pref_noad1 is same as ori_pref_noad
    #     df_tidy.loc[(df_tidy.cell_id == icell), 'ori_pref_noad2'] = ori_pref_run2[icell, 0] # take from another run

    #     df_tidy.loc[df_tidy['cell_id']==icell, 'filter_cell_well_fit_ad_250'] = well_fit_ad_250[icell]
    #     df_tidy.loc[df_tidy['cell_id']==icell, 'filter_cell_well_fit_ad_750'] = well_fit_ad_750[icell]
    
    # break

    ## save df_tidy as csv
    df_tidy.to_csv(os.path.join(dir_data, 'df_tidy_continuous_vis_pval.csv'), index=False)

# clear_output()


# # batch load df_tidy

# In[7]:


meta = pd.concat([meta_V1, meta_LM, meta_LI], axis=0).reset_index(drop=True)
# meta = meta_LI.copy()
# meta.sample(5, random_state=0)
meta


# In[8]:


nset = len(meta)
df_tidy = pd.DataFrame()

for iset in range(nset):
    print(f'iset={iset}, nset={nset}')

    mouse = meta.loc[iset, 'mouse'].astype(str)
    imouse = 'i' + mouse
    date = meta.loc[iset, 'date'].astype(str)
    area = meta.loc[iset, 'area']
    # sess = '00' + meta.loc[iset, 'num'].astype(str)
    sess = meta.loc[iset, 'num']
    print(imouse, date, area, sess)
    dir_identifier = f'{area}_{imouse}_{date}_{sess}'

    dir_data = r'\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images\mat_inter/'.replace('\\', '/')
    csv_filename = 'df_tidy_continuous_vis_pval.csv'
    try:
        df_tidy_now = pd.read_csv(os.path.join(dir_data, dir_identifier, csv_filename))
    except:
        dir_identifier = dir_identifier + '_cellpose'
        df_tidy_now = pd.read_csv(os.path.join(dir_data, dir_identifier, csv_filename))

    df_tidy_now['mouse'] = mouse
    df_tidy_now['date'] = date
    df_tidy_now['sess'] = sess
    df_tidy_now['cell_id'] = (df_tidy_now.date.astype(str) + '_' 
                               + df_tidy_now.sess.astype(str) + '_' 
                               + df_tidy_now.cell_id.astype(str)) # cell_id adjusted to be unique across mice, dates, sessions
    df_tidy = pd.concat([df_tidy, df_tidy_now], axis=0).reset_index(drop=True)
    # break
    
clear_output()


# In[9]:


print(df_tidy.mouse.unique(), 
      df_tidy.date.unique(), 
      df_tidy.sess.unique(), 
      df_tidy.area.unique(), 
      df_tidy.isi.unique(), 
      df_tidy.stim1_contrast.unique(), 
      df_tidy.stim2_id.unique(), 
      df_tidy.resp_id.unique())
df_tidy.sample(5, random_state=0)


# # adaptation by area

# In[120]:


## suppress FutureWarning
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

def adp_in_area(df_filter):

    gb_R1 = df_filter[(df_filter.stim2_id == 0) & (df_filter.resp_id == 'R1')].groupby(['area', 'cell_id']).sum()['dfof'].values
    gb_R2 = df_filter[(df_filter.stim2_id == 0) & (df_filter.resp_id == 'R2')].groupby(['area', 'cell_id']).sum()['dfof'].values
    gb_adp = (gb_R2 - gb_R1) / (gb_R1 + 1e-7)
    gb_adp = gb_adp[np.abs(gb_adp) < 2]

    return gb_adp

df_filter = df_tidy[(df_tidy.isi == 250)
                    & (df_tidy.filter_cell_vis == True)
                    ]
gb_adp_V1 = adp_in_area(df_filter[df_filter.area == 'V1'])
gb_adp_LM = adp_in_area(df_filter[df_filter.area == 'LM'])
gb_adp_LI = adp_in_area(df_filter[df_filter.area == 'LI'])

# plt.hist(gb_adp_V1, bins=50, alpha=.5, label='V1')
# plt.hist(gb_adp_LM, bins=50, alpha=.5, label='LM')
# plt.hist(gb_adp_LI, bins=50, alpha=.5, label='LI')
# plt.legend();

# adp_mean_arr = [np.mean(gb_adp_V1), np.mean(gb_adp_LM), np.mean(gb_adp_LI)]
# adp_med_arr = [np.median(gb_adp_V1), np.median(gb_adp_LM), np.median(gb_adp_LI)]
# adp_sem_arr = [np.std(gb_adp_V1) / np.sqrt(len(gb_adp_V1)), 
#                np.std(gb_adp_LM) / np.sqrt(len(gb_adp_LM)),
#                 np.std(gb_adp_LI) / np.sqrt(len(gb_adp_LI))]
# plt.errorbar([1, 2, 3], adp_mean_arr, yerr=adp_sem_arr, label='mean', alpha=.5)
# plt.errorbar([1, 2, 3], adp_med_arr, yerr=adp_sem_arr, label='median', alpha=.5)
# plt.xticks([1, 2, 3], ['V1', 'LM', 'LI']);
# plt.ylim([-1, 0]);
# plt.legend(frameon=False);

# sns.boxplot(data=[gb_adp_V1, gb_adp_LM, gb_adp_LI], palette='Set2', notch=True, showfliers=False) # dont show the outliers beyond the caps
ax = sns.violinplot(data=[gb_adp_V1, gb_adp_LM, gb_adp_LI], 
                    palette=['blue', 'orange', 'green'], 
                    )
plt.setp(ax.collections, alpha=.6)
plt.setp(ax.lines, alpha=.5);

## ncell in each area
xpos = -0.2
ypos = -2.8
plt.text(xpos, ypos, f'n={len(gb_adp_V1)}', fontsize=14)
plt.text(xpos+1, ypos, f'n={len(gb_adp_LM)}', fontsize=14)
plt.text(xpos+2, ypos, f'n={len(gb_adp_LI)}', fontsize=14)

## p value
_, p_kruskal = stats.kruskal(gb_adp_V1, gb_adp_LM, gb_adp_LI)
plt.text(xpos+1.8, 2.5, 'p={:.1e}'.format(p_kruskal), fontsize=14)

plt.xticks([0, 1, 2], ['V1', 'LM', 'LI'], fontsize=16);
plt.ylim([-3, 3]);
plt.ylabel('Adaptation index', fontsize=16);
sns.despine();
plt.tight_layout();

dir_fig = repo_dir + r'\results\tuning bias'.replace('\\', '/')
plt.savefig(os.path.join(dir_fig, 'adp_by_area_grat8_vis.pdf'))


# In[10]:


df_filter.groupby('area').cell_id.nunique()


# In[12]:


## kruskal wallis test
print(stats.kruskal(gb_adp_V1, gb_adp_LM, gb_adp_LI))

## mann whitney u test (Wilcoxon rank-sum test)
stats.mannwhitneyu(gb_adp_V1, gb_adp_LM), stats.mannwhitneyu(gb_adp_V1, gb_adp_LI), stats.mannwhitneyu(gb_adp_LM, gb_adp_LI)


# ## population level adp across area

# In[168]:


df_tidy['date_sess'] = df_tidy.date.astype(str) + '_' + df_tidy.sess.astype(str)
df_tidy.date_sess.unique()


# In[179]:


# ## suppress FutureWarning
# import warnings
# warnings.simplefilter(action='ignore', category=FutureWarning)

def adp_in_area_pop(df_filter):

    gb_R1 = df_filter[(df_filter.stim2_id == 0) & (df_filter.resp_id == 'R1')].groupby(['mouse']).sum()['dfof'].values
    gb_R2 = df_filter[(df_filter.stim2_id == 0) & (df_filter.resp_id == 'R2')].groupby(['mouse']).sum()['dfof'].values
    gb_adp = (gb_R2 - gb_R1) / (gb_R1 + 1e-7)
    # gb_adp = gb_adp[np.abs(gb_adp) < 2]

    return gb_adp

df_filter = df_tidy[(df_tidy.isi == 250)
                    & (df_tidy.filter_cell_vis == True)
                    ]
gb_adp_V1 = adp_in_area_pop(df_filter[df_filter.area == 'V1'])
gb_adp_LM = adp_in_area_pop(df_filter[df_filter.area == 'LM'])
gb_adp_LI = adp_in_area_pop(df_filter[df_filter.area == 'LI'])

ax = sns.violinplot(data=[gb_adp_V1, gb_adp_LM, gb_adp_LI], 
                    palette=['blue', 'orange', 'green'], 
                    )
plt.setp(ax.collections, alpha=.6)
plt.setp(ax.lines, alpha=.5);

## ncell in each area
xpos = -0.2
ypos = -2.8
plt.text(xpos, ypos, f'n={len(gb_adp_V1)}', fontsize=14)
plt.text(xpos+1, ypos, f'n={len(gb_adp_LM)}', fontsize=14)
plt.text(xpos+2, ypos, f'n={len(gb_adp_LI)}', fontsize=14)

## p value
_, p_kruskal = stats.kruskal(gb_adp_V1, gb_adp_LM, gb_adp_LI)
plt.text(xpos+1.8, 2.5, 'p={:.3f}'.format(p_kruskal), fontsize=14)

plt.xticks([0, 1, 2], ['V1', 'LM', 'LI'], fontsize=16);
# plt.ylim([-3, 3]);
plt.ylabel('Adaptation index', fontsize=16);
sns.despine();
plt.tight_layout();

dir_fig = repo_dir + r'\results\tuning bias'.replace('\\', '/')
# plt.savefig(os.path.join(dir_fig, 'adp_by_area_grat8_vis.pdf'))


# In[180]:


## kruskal wallis test
print(stats.kruskal(gb_adp_V1, gb_adp_LM, gb_adp_LI))

## mann whitney u test (Wilcoxon rank-sum test)
stats.mannwhitneyu(gb_adp_V1, gb_adp_LM), stats.mannwhitneyu(gb_adp_V1, gb_adp_LI), stats.mannwhitneyu(gb_adp_LM, gb_adp_LI)


# ## old version adaptation

# In[122]:


# ## suppress FutureWarning
# import warnings
# warnings.simplefilter(action='ignore', category=FutureWarning)

def old_adp_in_area(df_filter, area):

    gb_R1 = (df_filter[(df_filter.area == area) 
                & (df_filter.stim2_id == 0) 
                & (df_filter.resp_id == 'R1')
                ]
            .groupby(['trial_id', 'cell_id']) # each trial_id is a pair of R1-R2
            .dfof.apply(lambda x: x.values[0]) # take dfof value from the only row in each group
            )
    gb_R2 = (df_filter[(df_filter.area == area) 
                & (df_filter.stim2_id == 0) 
                & (df_filter.resp_id == 'R2')
                ]
            .groupby(['trial_id', 'cell_id'])
            .dfof.apply(lambda x: x.values[0])
            )

    gb_adp = (gb_R2 - gb_R1) / (gb_R1 + 1e-7)

    gb = (df_filter[(df_filter.area == area) 
                & (df_filter.stim2_id == 0) 
                & (df_filter.resp_id == 'R1')
                ]
            .groupby(['trial_id', 'cell_id'])
            .dfof.sum().to_frame() # not used, just providing index
            )
    gb['adp'] = gb_adp
    gb = gb.drop(columns='dfof').reset_index()
    adp_cell = gb.groupby('cell_id').adp.mean().values
    adp_cell = adp_cell[np.abs(adp_cell) < 2]

    return adp_cell


df_filter = df_tidy[(df_tidy.isi == 250)
                    & (df_tidy.filter_cell_vis == True)
                    ]
gb_adp_V1 = old_adp_in_area(df_filter, 'V1')
gb_adp_LM = old_adp_in_area(df_filter, 'LM')
gb_adp_LI = old_adp_in_area(df_filter, 'LI')


ax = sns.violinplot(data=[gb_adp_V1, gb_adp_LM, gb_adp_LI], 
                    palette=['blue', 'orange', 'green'], 
                    )
plt.setp(ax.collections, alpha=.6)
plt.setp(ax.lines, alpha=.5);

## ncell in each area
xpos = -0.2
ypos = -2.9
plt.text(xpos, ypos, f'n={len(gb_adp_V1)}', fontsize=14)
plt.text(xpos+1, ypos, f'n={len(gb_adp_LM)}', fontsize=14)
plt.text(xpos+2, ypos, f'n={len(gb_adp_LI)}', fontsize=14)

# ## p value
# _, p_kruskal = stats.kruskal(gb_adp_V1, gb_adp_LM, gb_adp_LI)
# plt.text(xpos+1.8, 2.5, 'p={:.1e}'.format(p_kruskal), fontsize=14)

plt.xticks([0, 1, 2], ['V1', 'LM', 'LI'], fontsize=16);
plt.ylim([-3, 3]);
plt.ylabel('Adaptation index', fontsize=16);
sns.despine();
plt.tight_layout();

dir_fig = repo_dir + r'\results\tuning bias'.replace('\\', '/')
plt.savefig(os.path.join(dir_fig, 'old_ver_adp_by_area_grat8_vis.pdf'))


# # tuning bias crude max-ori

# ## filter cell well-max

# In[13]:


df_nrep = (df_tidy[df_tidy.resp_id == 'R2'] # only R2 has diff ori
            [['dfof', 'cell_id', 'resp_id', 'isi', 'stim2_id']]
              .groupby(['cell_id', 'isi', 'stim2_id']).count() # count trials per cell, isi, ori
            )
plt.hist(df_nrep.values.flatten(), bins=20);
min(df_nrep.values.flatten())


# In[15]:


## construct tuning_vec column

max_ori_bootstrap = np.array([])
nboot = 50

for iboot in tqdm(range(nboot)):
    max_ori = (df_tidy[df_tidy.resp_id == 'R2'] # only R2 has diff ori
                [['dfof', 'cell_id', 'resp_id', 'isi', 'stim2_id']]
                .groupby(['cell_id', 'isi', 'stim2_id'])
                .sample(frac=.7, random_state=iboot) # sample x trials per group. min rep = 48
                .groupby(['cell_id', 'isi', 'stim2_id'])
                .agg({'dfof': 'mean'}) # aggregate resp by cell, isi, ori
                .groupby(['cell_id', 'isi']).apply(lambda x: np.argmax(x)) # get ori with max resp
                .reset_index() # NOTE: due to prev groupby, cell_id and isi are in order
                .rename(columns={0: 'max_ori'}).max_ori.values
                )
    max_ori_bootstrap = np.append(max_ori_bootstrap, max_ori)


# In[16]:


max_ori_bootstrap = max_ori_bootstrap.reshape(nboot, -1) # reshape max_ori_bootstrap to: nboot x (ncell x nisi)

## sort max_ori_bootstrap in each col (across boots)
max_ori_bootstrap = np.sort(max_ori_bootstrap, axis=0)
# sns.heatmap(max_ori_bootstrap[:, 10:30], cmap='viridis', cbar=True, annot=True);

## count the most freq value in each col (across boots)
max_ori_mode = stats.mode(max_ori_bootstrap, axis=0)[0][0]
max_ori_mode_freq = stats.mode(max_ori_bootstrap, axis=0)[1][0]
filter_cell_well_max = (max_ori_mode_freq >= nboot * 0.7) # NOTE: strictness of well-max can be adjusted here

# ## query values at ith and (100-i)th percentile to see if their diff > 1
# percentile = 5 # NOTE: strictness of well-max can be adjusted here
# max_ori_boot_low = np.percentile(max_ori_bootstrap, percentile, axis=0)
# max_ori_boot_high = np.percentile(max_ori_bootstrap, 100-percentile, axis=0)
# max_ori_boot_var = max_ori_boot_high - max_ori_boot_low
# max_ori_boot_var.shape # df_tidy.cell_id.nunique() x df_tidy.isi.nunique() = 160 x 3

# ## use max_ori_boot_var to get filter_cell_well_max
# filter_cell_well_max = ((max_ori_boot_var <= 1) | (max_ori_boot_var == 7)) # 7 is equivalent to 1, due to circularity. NOTE: strictness of well-max can be adjusted here too
# # plt.plot(max_ori_boot_var[10:30], color='k');
# # plt.plot(filter_cell_well_max[10:30], color='r');

## merge filter_cell_well_max with df
df_well_max = (df_tidy[df_tidy.resp_id == 'R2'][['cell_id', 'isi']]
                .groupby(['cell_id', 'isi'])
                .first() # get first trial per cell, isi
                .reset_index() # due to prev groupby, cell_id and isi are in the same order as max_ori
                )
df_well_max['filter_cell_well_max'] = filter_cell_well_max
# df_well_max.filter_cell_well_max.sum() / df_well_max.shape[0]
df_well_max


# In[17]:


## inherit df_well_max['filter_cell_well_max']
df_tidy = df_tidy.merge(
      df_well_max, 
      on=['cell_id', 'isi'], how='left') # take filter_cell_well_max from df_well_max, for each cell and isi
df_tidy


# ## pref_ori for cell, no adapter or with adapter
# crude preference: take max-resp orientation as the preferred orien

# In[19]:


pref_ori_noad = (df_tidy[(df_tidy.resp_id == 'R2') & (df_tidy.isi == 6000)]
                .groupby(['cell_id', 'stim2_id'])
                [['dfof']].sum().reset_index()
                .groupby('cell_id') # NOTE: stim2_id is sorted due to prev groupby
                ['dfof'].apply(list).to_frame()
                .dfof.apply(lambda x: np.argmax(x)).to_frame()
                .rename(columns={'dfof': 'pref_ori_noad'})
                )

pref_ori_ad_250 = (df_tidy[(df_tidy.resp_id == 'R2') & (df_tidy.isi == 250)]
                .groupby(['cell_id', 'stim2_id'])
                [['dfof']].sum().reset_index()
                .groupby('cell_id') # NOTE: stim2_id is sorted due to prev groupby
                ['dfof'].apply(list).to_frame()
                .dfof.apply(lambda x: np.argmax(x)).to_frame()
                .rename(columns={'dfof': 'pref_ori_ad_250'})
                )

pref_ori_ad_750 = (df_tidy[(df_tidy.resp_id == 'R2') & (df_tidy.isi == 750)]
                .groupby(['cell_id', 'stim2_id'])
                [['dfof']].sum().reset_index()
                .groupby('cell_id') # NOTE: stim2_id is sorted due to prev groupby
                ['dfof'].apply(list).to_frame()
                .dfof.apply(lambda x: np.argmax(x)).to_frame()
                .rename(columns={'dfof': 'pref_ori_ad_750'})
                )

df_tidy = (df_tidy.merge(pref_ori_noad, on='cell_id')
            .merge(pref_ori_ad_250, on='cell_id')
            .merge(pref_ori_ad_750, on='cell_id'))
df_tidy


# In[20]:


def bin_ori(x):
    # bin pref ori to 0, 45, 90
    if x < 30:
        return 0
    elif x <= 60:
        return 45
    else:
        return 90

def distace_from_adapter(x):
    # adapter is always 0 deg
    x = 22.5*x # convert from stim id int to degree
    if x > 90:
        x = 180 - x
    return x


df_tidy['pref_unadapted_distance'] = df_tidy['pref_ori_noad'].apply(lambda x: distace_from_adapter(x)) # unadapted pref ori, distance from adapter
df_tidy['pref_unadapted_distance_bin'] = df_tidy['pref_unadapted_distance'].apply(lambda x: bin_ori(x)) # bin the distance to 0, 45, 90

df_tidy.loc[(df_tidy.isi == 250), 'pref_adapted_distance'] = df_tidy['pref_ori_ad_250'].apply(lambda x: distace_from_adapter(x)) # isi 250 adapted pref ori, distance from adapter # TODO: refactor pref_ori_ad_250 similarly, so it's only one column of pref_ori_ad, but can filter by isi
df_tidy.loc[(df_tidy.isi == 750), 'pref_adapted_distance'] = df_tidy['pref_ori_ad_750'].apply(lambda x: distace_from_adapter(x)) # isi 750 adapted pref ori, distance from adapter

df_tidy['tuning_bias'] = df_tidy['pref_adapted_distance'] - df_tidy['pref_unadapted_distance'] # distance from adapter, adapted - unadapted. if tuning repelled from adapter, this is positive; attracted, negative

# dir_data = r'\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images\mat_inter\df tidy for plots/'.replace('\\', '/')
# df_tidy.to_csv(os.path.join(dir_data, 'df_tidy_tuning_bias_nofit.csv'), index=False)

df_tidy.sample(5)


# ## tuning bias plot no fit

# In[124]:


df_filter = df_tidy[(df_tidy.filter_cell_vis == True) 
                    & (df_tidy.filter_cell_well_max == True)
                    # & (df_tidy.filter_cell_stim == True)
                    & (df_tidy.isi == 250)
                    ]
df_filter = df_filter.groupby('cell_id').first().reset_index() # drop duplicate cell_id
df_filter.groupby('area').cell_id.nunique()


# In[151]:


fig, axes = plt.subplots(1, 1, figsize=(6, 5), sharey=True)
ax = sns.pointplot(x='pref_unadapted_distance_bin', y='tuning_bias', data=df_filter, hue='area',
                   dodge=True, markers='.', 
                   errorbar=('ci', 68), #errwidthfloat=1, capsize=.1,
                   )

ax.text(0.15, 0.99, '***', transform=ax.transAxes, fontsize=14)
ax.text(0.5, 0.99, 'n.s.', transform=ax.transAxes, fontsize=14)
ax.text(0.8, 0.99, 'n.s.', transform=ax.transAxes, fontsize=14)

## ncell in each area for each bin
jitter_large = 0.3
jitter_small = 0.1
y_text = -14
ax.text(-jitter_large, y_text, f'V1 n={df_filter[(df_filter.area == "V1") & (df_filter.pref_unadapted_distance_bin == 0)].cell_id.nunique()}', fontsize=12)
ax.text(1-jitter_small, y_text, f'{df_filter[(df_filter.area == "V1") & (df_filter.pref_unadapted_distance_bin == 45)].cell_id.nunique()}', fontsize=12)
ax.text(2-jitter_small, y_text, f'{df_filter[(df_filter.area == "V1") & (df_filter.pref_unadapted_distance_bin == 90)].cell_id.nunique()}', fontsize=12)
ax.text(-jitter_large, y_text-2, f'LM n={df_filter[(df_filter.area == "LM") & (df_filter.pref_unadapted_distance_bin == 0)].cell_id.nunique()}', fontsize=12)
ax.text(1-jitter_small, y_text-2, f'{df_filter[(df_filter.area == "LM") & (df_filter.pref_unadapted_distance_bin == 45)].cell_id.nunique()}', fontsize=12)
ax.text(2-jitter_small, y_text-2, f'{df_filter[(df_filter.area == "LM") & (df_filter.pref_unadapted_distance_bin == 90)].cell_id.nunique()}', fontsize=12)
ax.text(-jitter_large, y_text-4, f'LI n={df_filter[(df_filter.area == "LI") & (df_filter.pref_unadapted_distance_bin == 0)].cell_id.nunique()}', fontsize=12)
ax.text(1-jitter_small, y_text-4, f'{df_filter[(df_filter.area == "LI") & (df_filter.pref_unadapted_distance_bin == 45)].cell_id.nunique()}', fontsize=12)
ax.text(2-jitter_small, y_text-4, f'{df_filter[(df_filter.area == "LI") & (df_filter.pref_unadapted_distance_bin == 90)].cell_id.nunique()}', fontsize=12)

plt.setp(ax.collections, alpha=.5);
plt.setp(ax.lines, alpha=.5);
plt.axhline(0, color='gray', linestyle='-', alpha=.3);
plt.legend(frameon=False, loc='right');
plt.xlabel('Preferrence distance from adapter', fontsize=16);
plt.ylabel('Tuning bias', fontsize=16);
plt.ylim([-20, 32]);
sns.despine();
plt.tight_layout();

plt.savefig(os.path.join(dir_fig, 'tuning_bias_nofit_by_area_grat8_vis_wellmax_ncell.pdf'))


# In[133]:


(df_filter
    # .sort_values(by='area')
    # .sort_values(by='pref_unadapted_distance_bin')
    .groupby(['area', 'pref_unadapted_distance_bin'], 
            #  sort=False
             )
    .cell_id.nunique())


# In[53]:


bias_0_V1 = df_filter[(df_filter.pref_unadapted_distance_bin == 0) & (df_filter.area == 'V1')].tuning_bias.values
bias_0_LM = df_filter[(df_filter.pref_unadapted_distance_bin == 0) & (df_filter.area == 'LM')].tuning_bias.values
bias_0_LI = df_filter[(df_filter.pref_unadapted_distance_bin == 0) & (df_filter.area == 'LI')].tuning_bias.values
print(stats.kruskal(bias_0_V1, bias_0_LM, bias_0_LI))

bias_45_V1 = df_filter[(df_filter.pref_unadapted_distance_bin == 45) & (df_filter.area == 'V1')].tuning_bias.values
bias_45_LM = df_filter[(df_filter.pref_unadapted_distance_bin == 45) & (df_filter.area == 'LM')].tuning_bias.values
bias_45_LI = df_filter[(df_filter.pref_unadapted_distance_bin == 45) & (df_filter.area == 'LI')].tuning_bias.values
print(stats.kruskal(bias_45_V1, bias_45_LM, bias_45_LI))

bias_90_V1 = df_filter[(df_filter.pref_unadapted_distance_bin == 90) & (df_filter.area == 'V1')].tuning_bias.values
bias_90_LM = df_filter[(df_filter.pref_unadapted_distance_bin == 90) & (df_filter.area == 'LM')].tuning_bias.values
bias_90_LI = df_filter[(df_filter.pref_unadapted_distance_bin == 90) & (df_filter.area == 'LI')].tuning_bias.values
print(stats.kruskal(bias_90_V1, bias_90_LM, bias_90_LI))


# In[23]:


fig, axes = plt.subplots(1, 1, figsize=(6, 5), sharey=True)
ax = sns.pointplot(x='pref_unadapted_distance', y='tuning_bias', data=df_filter, hue='area',
                   dodge=True, markers='.', 
                   errorbar=('ci', 68), #errwidthfloat=1, capsize=.1,
                   )
plt.setp(ax.collections, alpha=.5);
plt.setp(ax.lines, alpha=.5);
plt.axhline(0, color='gray', linestyle='-', alpha=.3);
plt.legend(frameon=False);


# In[26]:


df_filter.groupby(['area', 'pref_unadapted_distance']).cell_id.nunique()


# In[ ]:


# df_filter.groupby('pref_unadapted_distance').cell_id.nunique().sort_index().plot(kind='bar')
# # plt.xticks(np.arange(0, 5, 1), np.arange(0, 90+22.5, 22.5), rotation=45)
# plt.xlabel('pref ori (no adapter) degree')
# plt.ylabel('number of cells');


# # tuning bias preprocessing
# x: distance(pref_ori_unadapted, 0 deg adapter_ori).binned  
# y: distance(pref_ori_unadapted, 0 deg) - distance(pref_ori_adapted, 0 deg)

# ## pref_ori for cell & isi
# ~~crude preference: take max-resp orientation as the preferred orien~~  
# fitted preference: use pref from von mises curve fit

# In[10]:


def bin_ori(x):
    # bin pref ori to 0, 45, 90
    if x < 30:
        return 0
    elif x <= 60:
        return 45
    else:
        return 90

def bin_ori_finer(x):
    # bin pref ori to n bins, with equal bin width
    nbin = 4 # n_edge = nbin + 1
    bin_width = 90 / nbin
    return (x // bin_width) * bin_width

def distance_from_adapter(x):
    # adapter is always 0 deg
    # x = 22.5*x # convert from stim id int to degree
    if x > 90:
        x = 180 - x
    return x


df_tidy['pref_unadapted_distance'] = df_tidy['ori_pref_noad'].apply(lambda x: distance_from_adapter(x)) # unadapted pref ori, distance from adapter
df_tidy['pref_unadapted_distance_bin'] = df_tidy['pref_unadapted_distance'].apply(lambda x: bin_ori(x)) # bin the distance to 0, 45, 90
# df_tidy['pref_unadapted_distance_bin'] = df_tidy['pref_unadapted_distance'].apply(lambda x: bin_ori_finer(x)) # bin the distance to 5 edges, 0, 22.5, 45, 67.5, 90

df_tidy.loc[(df_tidy.isi == 250), 'pref_adapted_distance'] = df_tidy['ori_pref_ad_250'].apply(lambda x: distance_from_adapter(x)) # isi 250 adapted pref ori, distance from adapter # TODO: refactor pref_ori_ad_250 similarly, so it's only one column of pref_ori_ad, but can filter by isi
df_tidy.loc[(df_tidy.isi == 750), 'pref_adapted_distance'] = df_tidy['ori_pref_ad_750'].apply(lambda x: distance_from_adapter(x)) # isi 750 adapted pref ori, distance from adapter
df_tidy['tuning_bias'] = df_tidy['pref_adapted_distance'] - df_tidy['pref_unadapted_distance'] # distance from adapter, adapted - unadapted. if tuning repelled from adapter, this is positive; attracted, negative

df_tidy.sample(5, random_state=0)


# In[11]:


df_control = df_tidy.copy()

df_control['pref_unadapted_distance'] = df_control['ori_pref_noad1'].apply(lambda x: distance_from_adapter(x)) # unadapted pref ori, distance from adapter
df_control['pref_unadapted_distance_bin'] = df_control['pref_unadapted_distance'].apply(lambda x: bin_ori(x)) # bin the distance to 0, 45, 90
# df_control['pref_unadapted_distance_bin'] = df_control['pref_unadapted_distance'].apply(lambda x: bin_ori_finer(x))

df_control['pref_adapted_distance'] = df_control['ori_pref_noad2'].apply(lambda x: distance_from_adapter(x))
df_control['tuning_bias'] = df_control['pref_adapted_distance'] - df_control['pref_unadapted_distance']

df_control.sample(5, random_state=0)


# ## merge df real vs control

# In[12]:


assert len(df_tidy) == len(df_control) # same number of rows
assert (df_tidy.columns == df_control.columns).all() # all columns are the same
assert np.sum(df_control['pref_unadapted_distance_bin'].values == df_tidy['pref_unadapted_distance_bin'].values) == len(df_tidy) # same pref_unadapted_distance_bin due to ori_pref_noad1 == ori_pref_noad
assert np.sum(df_tidy.tuning_bias.values == df_control.tuning_bias.values) < len(df_tidy) # different tuning bias
print('for tuning bias plot, only tuning bias col differs between control and real df')


# In[13]:


df_tidy['tuning_bias_control'] = df_control['tuning_bias']
df_tidy


# ## save & reload

# In[15]:


dir_df = r'\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images\mat_inter\df tidy for plots/'.replace('\\', '/')
csv_filename = dir_df + 'df_tidy_tuning_bias_diff_exclude_bad_LI.csv'

# chunks = np.array_split(df_tidy.index, 100) # split into 100 chunks
# for chunk, subset in enumerate(tqdm(chunks)):
#     if chunk == 0: # first row
#         df_tidy.loc[subset].to_csv(csv_filename, mode='w', index=True)
#     else:
#         df_tidy.loc[subset].to_csv(csv_filename, header=None, mode='a', index=True)

df_tidy = pd.read_csv(os.path.join(dir_df, csv_filename), index_col=0) # this csv excludes bad LI data
df_tidy


# ## if use all trials, not half trials

# In[31]:


df_tidy_all_trials = pd.read_csv(os.path.join(dir_df, 'df_tidy_tuning_bias.csv'), index_col=0)
df_tidy_control_all_trials = pd.read_csv(os.path.join(dir_df, 'df_tidy_tuning_bias_control.csv'), index_col=0)


# In[36]:


assert len(df_tidy_all_trials) == len(df_tidy_control_all_trials) # same number of rows
assert (df_tidy_all_trials.columns == df_tidy_control_all_trials.columns).all() # all columns are the same
# assert np.sum(df_tidy_control_all_trials['pref_unadapted_distance_bin'].values == df_tidy_all_trials['pref_unadapted_distance_bin'].values) == len(df_tidy_all_trials) # same pref_unadapted_distance_bin due to ori_pref_noad1 == ori_pref_noad
assert np.sum(df_tidy_all_trials.tuning_bias.values == df_tidy_control_all_trials.tuning_bias.values) < len(df_tidy_all_trials) # different tuning bias
print('for tuning bias plot, only tuning bias and pref_unadapted_distance_bin differs between control and real df')
df_tidy_all_trials['tuning_bias_control'] = df_tidy_control_all_trials['tuning_bias']
df_tidy_all_trials['pref_unadapted_distance_bin_control'] = df_tidy_control_all_trials['pref_unadapted_distance_bin']
df_tidy_all_trials


# In[40]:


df_filter = df_tidy_all_trials[(df_tidy_all_trials.trial_id > -1) # placeholder, always true
                            & (df_tidy_all_trials.filter_cell_vis == True)
                            & (df_tidy_all_trials.filter_cell_well_fit == True)
                            & (df_tidy_all_trials.isi == 250)
                            ]
df_filter = df_filter.groupby('cell_id')[['area', 'pref_unadapted_distance_bin', 'pref_unadapted_distance_bin_control', 'tuning_bias', 'tuning_bias_control']].first().reset_index()
df_filter['tuning_bias_diff'] = df_filter['tuning_bias'] - df_filter['tuning_bias_control']
df_filter


# ### tuning bias plot (all trials)

# In[41]:


fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
jitter = 1
for i, iarea in enumerate(df_filter.area.unique()):
    axes[0].errorbar(x=np.array(sorted(df_filter[df_filter.area == iarea].pref_unadapted_distance_bin.unique())) + i * jitter,
                    y=df_filter[df_filter.area == iarea].groupby('pref_unadapted_distance_bin').tuning_bias.mean(), 
                    yerr=df_filter[df_filter.area == iarea].groupby('pref_unadapted_distance_bin').tuning_bias.sem(), 
                    capsize=5, alpha=.7)


    axes[1].errorbar(x=np.array(sorted(df_filter[df_filter.area == iarea].pref_unadapted_distance_bin_control.unique())) + i * jitter,
                    y=df_filter[df_filter.area == iarea].groupby('pref_unadapted_distance_bin').tuning_bias_control.mean(),
                    yerr=df_filter[df_filter.area == iarea].groupby('pref_unadapted_distance_bin').tuning_bias_control.sem(),
                    capsize=5, alpha=.7)
    
    axes[2].errorbar(x=np.array(sorted(df_filter[df_filter.area == iarea].pref_unadapted_distance_bin.unique())) + i * jitter,
                    y=df_filter[df_filter.area == iarea].groupby('pref_unadapted_distance_bin').tuning_bias.mean() - df_filter[df_filter.area == iarea].groupby('pref_unadapted_distance_bin').tuning_bias_control.mean(),
                    yerr=df_filter[df_filter.area == iarea].groupby('pref_unadapted_distance_bin').tuning_bias_diff.sem(),
                    capsize=5, alpha=.7)

axes[0].axhline(0, color='gray', linestyle='-', alpha=.3);
axes[1].axhline(0, color='gray', linestyle='-', alpha=.3);
axes[2].axhline(0, color='gray', linestyle='-', alpha=.3);

axes[0].set_title('real')
axes[0].set_xlabel('pref_unadapted_distance_bin')
axes[0].set_ylabel('tuning_bias')
axes[0].set_xticks([0, 45, 90])
axes[0].set_xticklabels([0, 45, 90])

axes[1].set_title('control');
axes[1].set_xticks([0, 45, 90])
axes[1].set_xticklabels([0, 45, 90])

axes[2].set_title('diff');
axes[2].set_xticks([0, 45, 90])
axes[2].set_xticklabels([0, 45, 90]);


# ### stats sig

# In[50]:


pref_unadapted_distance_bin = 90
data = np.array([
                df_filter[(df_filter.area == 'V1') & (df_filter.pref_unadapted_distance_bin == pref_unadapted_distance_bin)].tuning_bias_diff.values,
                df_filter[(df_filter.area == 'LM') & (df_filter.pref_unadapted_distance_bin == pref_unadapted_distance_bin)].tuning_bias_diff.values,
                df_filter[(df_filter.area == 'LI') & (df_filter.pref_unadapted_distance_bin == pref_unadapted_distance_bin)].tuning_bias_diff.values,
                ], dtype=object)
data.shape, data[0].shape, data[1].shape, data[2].shape


# In[51]:


## t test
from scipy.stats import ttest_ind

## assumption: equal variance
print(np.var(data[0]), np.var(data[1]),) # assumption not met

t, p = ttest_ind(data[0], data[1], equal_var=False, alternative='less') #, permutations=10000, random_state=0)
t, p


# In[52]:


## one way ANOVA
from scipy.stats import f_oneway

## assumption: equal variance
print(np.var(data[0]), np.var(data[1]), np.var(data[2]))
## assumption: normality
from scipy.stats import shapiro
print(shapiro(data[0]).pvalue, shapiro(data[1]).pvalue, shapiro(data[2]).pvalue)

f_oneway(data[0], data[1], data[2])


# In[53]:


from scipy.stats import kruskal
kruskal(data[0], data[1], data[2])


# # tuning bias plot

# In[16]:


df_filter = df_tidy[(df_tidy.trial_id > -1) # placeholder, always true
                    & (df_tidy.filter_cell_vis == True)
                  #   & (df_tidy.filter_cell_stim == True)
                    & (df_tidy.filter_cell_well_fit == True)
                    # & (df_tidy.filter_cell_well_fit_ad_250 == True)
                    & (df_tidy.isi == 250)
                    # & (df_tidy.filter_cell_well_fit_ad_750 == True)
                    # & (df_tidy.isi == 750)
                    ]
df_filter = df_filter.groupby('cell_id')[['area', 'pref_unadapted_distance_bin', 'tuning_bias', 'tuning_bias_control']].first().reset_index()
df_filter['tuning_bias_diff'] = df_filter['tuning_bias'] - df_filter['tuning_bias_control']
df_filter


# In[17]:


# set plt style to default
plt.style.use('default')

# set font size
plt.rcParams.update({'font.size': 14})

fig, axes = plt.subplots(1, 3, figsize=(10, 4), sharey=True, sharex=True)
jitter = 1
for i, iarea in enumerate(df_filter.area.unique()):
    axes[0].errorbar(x=np.array(sorted(df_filter[df_filter.area == iarea].pref_unadapted_distance_bin.unique())) + i * jitter,
                    y=df_filter[df_filter.area == iarea].groupby('pref_unadapted_distance_bin').tuning_bias.mean(), 
                    yerr=df_filter[df_filter.area == iarea].groupby('pref_unadapted_distance_bin').tuning_bias.sem(), 
                    capsize=5, alpha=.7)


    axes[1].errorbar(x=np.array(sorted(df_filter[df_filter.area == iarea].pref_unadapted_distance_bin.unique())) + i * jitter,
                    y=df_filter[df_filter.area == iarea].groupby('pref_unadapted_distance_bin').tuning_bias_control.mean(),
                    yerr=df_filter[df_filter.area == iarea].groupby('pref_unadapted_distance_bin').tuning_bias_control.sem(),
                    capsize=5, alpha=.7)
    
    axes[2].errorbar(x=np.array(sorted(df_filter[df_filter.area == iarea].pref_unadapted_distance_bin.unique())) + i * jitter,
                    y=df_filter[df_filter.area == iarea].groupby('pref_unadapted_distance_bin').tuning_bias_diff.mean(),
                    yerr=df_filter[df_filter.area == iarea].groupby('pref_unadapted_distance_bin').tuning_bias_diff.sem(),
                    capsize=5, alpha=.7)

axes[0].axhline(0, color='gray', linestyle='-', alpha=.3);
axes[1].axhline(0, color='gray', linestyle='-', alpha=.3);
axes[2].axhline(0, color='gray', linestyle='-', alpha=.3);

axes[0].set_title('real')
axes[0].set_xlabel('pref_unadapted_distance_bin')
axes[0].set_ylabel('tuning_bias')
axes[0].set_xticks([0, 45, 90])
axes[0].set_xticklabels([0, 45, 90])

axes[1].set_title('control');
axes[1].set_xticks([0, 45, 90])
axes[1].set_xticklabels([0, 45, 90])

axes[2].set_title('diff');
axes[2].set_xticks([0, 45, 90])
axes[2].set_xticklabels([0, 45, 90]);

# fig.savefig(os.path.join(r'C:\Users\ll357\Documents\inter\results\tuning bias'.replace('\\', '/'), 'tuning_bias_by_area_with_control.pdf'), bbox_inches='tight')


# In[18]:


fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
g1 = sns.pointplot(x='pref_unadapted_distance_bin', y='tuning_bias', 
                   data=df_filter[df_filter.area != 'LI'], hue='area',
                   errorbar="se", errwidthfloat=1, capsize=.1,
                   ax=axes[0],
                   )
g2 = sns.pointplot(x='pref_unadapted_distance_bin', y='tuning_bias_control', 
                   data=df_filter[df_filter.area != 'LI'], hue='area',
                   errorbar="se", errwidthfloat=1, capsize=.1,
                   ax=axes[1],
                  )
axes[0].axhline(0, color='gray', linestyle='-', alpha=.3);
axes[1].axhline(0, color='gray', linestyle='-', alpha=.3);

g2.set(ylabel=None) # remove ylabel
plt.setp(axes[0].collections, alpha=.5);
plt.setp(axes[1].collections, alpha=.5);
plt.setp(axes[0].lines, alpha=.7);
plt.setp(axes[1].lines, alpha=.7);
fig.tight_layout();


# In[19]:


fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
g1 = sns.pointplot(x='pref_unadapted_distance_bin', y='tuning_bias', data=df_filter, hue='area',
                   errorbar="se", errwidthfloat=1, capsize=.1,
                   ax=axes[0], dodge=True,
                   )
g2 = sns.pointplot(x='pref_unadapted_distance_bin', y='tuning_bias_control', data=df_filter, hue='area',
                  errorbar="se", errwidthfloat=1, capsize=.1,
                  ax=axes[1], dodge=True,
                  )
axes[0].axhline(0, color='gray', linestyle='-', alpha=.3);
axes[1].axhline(0, color='gray', linestyle='-', alpha=.3);

g2.set(ylabel=None) # remove ylabel
plt.setp(axes[0].collections, alpha=.5);
plt.setp(axes[1].collections, alpha=.5);
plt.setp(axes[0].lines, alpha=.7);
plt.setp(axes[1].lines, alpha=.7);
fig.tight_layout();

# plt.savefig(os.path.join(dir_fig, 'tuning_bias_by_area.pdf'))


# In[20]:


fig, axes = plt.subplots(1, 1, figsize=(9, 5), sharey=True)
g1 = sns.pointplot(x='pref_unadapted_distance_bin', y='tuning_bias_diff', data=df_filter, hue='area',
                   errorbar='se', 
                #    errorbar=('ci', 68), 
                   errwidthfloat=1, capsize=.1,
                   ax=axes, dodge=True,
                   )

# annotation above each dot, ncell
ncell_bin = df_filter.groupby('pref_unadapted_distance_bin').cell_id.nunique().sort_index().values
ylim = axes.get_ylim()
axes.set_ylim(ylim[0], ylim[1] + 1.5)
for i in range(len(ncell_bin)):
    axes.annotate(str(ncell_bin[i]), (i, ylim[1] + 0.1), ha='center', va='center', size=12)

axes.legend(frameon=False, bbox_to_anchor=(1.5, 0.5))
g1.set(xlabel='pref ori distance from 0 deg adapter')
g1.set(ylabel='tuning bias diff (real - control)')
axes.axhline(0, color='gray', linestyle='-', alpha=.3);
plt.setp(axes.collections, alpha=.5);
plt.setp(axes.lines, alpha=.7);
fig.tight_layout();

# plt.savefig(os.path.join(dir_fig, 'tuning_bias_control_diff_by_area.pdf'))


# ### stats sig

# In[25]:


pref_unadapted_distance_bin = 90
data = np.array([
                df_filter[(df_filter.area == 'V1') & (df_filter.pref_unadapted_distance_bin == pref_unadapted_distance_bin)].tuning_bias_diff.values,
                df_filter[(df_filter.area == 'LM') & (df_filter.pref_unadapted_distance_bin == pref_unadapted_distance_bin)].tuning_bias_diff.values,
                df_filter[(df_filter.area == 'LI') & (df_filter.pref_unadapted_distance_bin == pref_unadapted_distance_bin)].tuning_bias_diff.values,
                ], dtype=object)
data.shape, data[0].shape, data[1].shape, data[2].shape


# In[26]:


## t test
from scipy.stats import ttest_ind

## assumption: equal variance
print(np.var(data[0]), np.var(data[1]),) # assumption not met

t, p = ttest_ind(data[0], data[1], equal_var=False, alternative='less') #, permutations=10000, random_state=0)
t, p


# In[27]:


## one way ANOVA
from scipy.stats import f_oneway

## assumption: equal variance
print(np.var(data[0]), np.var(data[1]), np.var(data[2]))
## assumption: normality
from scipy.stats import shapiro
print(shapiro(data[0]).pvalue, shapiro(data[1]).pvalue, shapiro(data[2]).pvalue)

f_oneway(data[0], data[1], data[2])


# In[28]:


from scipy.stats import kruskal
kruskal(data[0], data[1], data[2])


# ### ncell by area & well fit %

# In[29]:


df_filter.groupby(['pref_unadapted_distance_bin', 'area']).cell_id.nunique(), \
df_tidy.groupby(['pref_unadapted_distance_bin', 'area']).cell_id.nunique(), 


# In[63]:


# stacked bar plot, colored by well_fit, for each area
df_tidy.groupby(['pref_unadapted_distance_bin', 'area', 'filter_cell_well_fit']).cell_id.nunique().unstack().plot(kind='bar', stacked=True, figsize=(12, 5), title='well_fit', colormap='tab10', legend='reverse', width=0.7);

# set legend frameon=False
plt.legend(frameon=False);

# set legend text content
plt.legend(['not well fit', 'well fit'], frameon=False);

plt.xlabel('preferred orientation distance from 0 degree adapter, area')
plt.ylabel('cell count')
plt.title('');
plt.tight_layout();

# plt.savefig(os.path.join(dir_fig, 'well_fit_by_area.pdf'))


# In[66]:


# stacked bar plot, colored by well_fit, for each area
(df_tidy.groupby(['pref_unadapted_distance_bin', 'area', 'filter_cell_well_max'])
 .cell_id.nunique().unstack()
 .plot(kind='bar', stacked=True, 
       figsize=(12, 5), title='well_max', 
       colormap='tab10', legend='reverse', width=0.7))

# set legend frameon=False
plt.legend(frameon=False);

# set legend text content
plt.legend(['not well max', 'well max'], frameon=False);

plt.xlabel('preferred orientation distance from 0 degree adapter, area')
plt.ylabel('cell count')
plt.title('');
plt.tight_layout();

plt.savefig(os.path.join(dir_fig, 'well_max_by_area.pdf'))


# # df for tuning curve

# In[154]:


df_tuning = df_tidy[['dfof', 'cell_id', 'resp_id', 'isi', 'stim2_id', 
                    'area', 'filter_cell_vis', 'filter_cell_well_fit', 'filter_cell_well_fit_ad_250', 
                    'ori_pref_ad_250', 'ori_pref_noad', 
                    'pref_unadapted_distance', 'pref_unadapted_distance_bin', 'pref_adapted_distance', 
                    'tuning_bias', #'tuning_bias_control',
                    ]]
df_tuning = df_tuning[(df_tuning.isi > -1) # placeholder, always true
                    # & (df_tuning.area == 'LI')
                    & (df_tuning.filter_cell_vis == True)
                    & (df_tuning.filter_cell_well_fit == True)
                    # & (df_tuning.filter_cell_well_fit_ad_250 == True) # only use cells that are well fit in both conditions: noad and ad 250
                    & ((df_tuning.isi == 250) | (df_tuning.isi == 6000))
                    # & (df_tuning.pref_unadapted_distance_bin == 90)
                    ]
                    
df_tuning['tuning_noad'] = np.pi
df_tuning['tuning_250'] = np.pi

for icell in tqdm(df_tuning.cell_id.unique()):
    ## tuning curve when isi = 6000, no adapter
    tuning_noad = df_tuning.loc[(df_tuning.cell_id == icell) & (df_tuning.isi == 6000) & (df_tuning.resp_id == 'R2'), :].groupby(['stim2_id'], sort=True).agg(np.nanmean)['dfof'].values# groupby sort: sorted by key, aka stim2_id. take R2 of no adapter trials to get tuning curve when no adapter
    df_tuning.loc[(df_tuning.cell_id == icell), 'tuning_noad'] = df_tuning.loc[(df_tuning.cell_id == icell), 'tuning_noad'].apply(lambda x: tuning_noad)
    
    ## tuning curve when isi = 250
    tuning_250 = df_tuning.loc[(df_tuning.cell_id == icell) & (df_tuning.isi == 250) & (df_tuning.resp_id == 'R2'), :].groupby(['stim2_id'], sort=True).agg(np.nanmean)['dfof'].values
    df_tuning.loc[(df_tuning.cell_id == icell), 'tuning_250'] = df_tuning.loc[(df_tuning.cell_id == icell), 'tuning_250'].apply(lambda x: tuning_250)

df_tuning.sort_values(by=['tuning_bias'], inplace=True) # df_tuning sort by tuning_bias_distance
df_tuning.sample(5, random_state=0) # NOTE: where isi=6000, pref_adapted_distance and tuning_bias are NaN


# ## polar plot of tuning
# before and after adaptation  
# filter cells 

# In[161]:


for icell in tqdm(df_tuning[df_tuning.area != 'V1'].cell_id.unique()):

    gOSI_noad = df_tuning[df_tuning.cell_id == icell].gOSI_noad.values[0]
    gOSI_250 = df_tuning[df_tuning.cell_id == icell].gOSI_250.values[0]

    tuning_noad = df_tuning[df_tuning.cell_id == icell].tuning_noad.values[0]
    tuning_noad = np.append(tuning_noad, tuning_noad) # repeat 8 values twice to make 16 values for polar plot
    tuning_noad = np.append(tuning_noad, tuning_noad[0]) # repeat first value at the end to close the circle
    ori_pref_noad = df_tuning[df_tuning.cell_id == icell].ori_pref_noad.values[0]
    # print('ori_pref_noad deg: ', np.round(ori_pref_noad, 2))
    ori_pref_noad = ori_pref_noad * np.pi / 180 # degree to radian

    tuning_250 = df_tuning[df_tuning.cell_id == icell].tuning_250.values[0]
    tuning_250 = np.append(tuning_250, tuning_250)
    tuning_250 = np.append(tuning_250, tuning_250[0])
    ori_pref_ad_250 = df_tuning[df_tuning.cell_id == icell].ori_pref_ad_250.values[0]
    # print('ori_pref_ad_250 deg: ', np.round(ori_pref_ad_250, 2))
    ori_pref_ad_250 = ori_pref_ad_250 * np.pi / 180

    tuning_bias = df_tuning[df_tuning.cell_id == icell].tuning_bias.values[0]
    # print('no adapter: ', np.round(tuning_noad, 2))
    # print('250 ms: ', np.round(tuning_250, 2))
    # print('tuning_bias: ', np.round(tuning_bias, 2))

    fig, ax = plt.subplots(1, 1, figsize=(12, 12), subplot_kw=dict(projection='polar'))
    ax.plot(np.linspace(0, 2*np.pi, 17), tuning_noad, alpha=.5, linewidth=10)
    ax.plot(np.linspace(0, 2*np.pi, 17), tuning_250, alpha=.5, linewidth=10)

    # plot adapter ori as a line
    min_val = np.min([np.min(tuning_noad), np.min(tuning_250)])
    max_val = np.max([np.max(tuning_noad), np.max(tuning_250)])
    # ax.plot([0, np.pi], [max_val, max_val], color='gray', linewidth=5, alpha=.5, label='adapter ori')

    # ## plot 90 deg as a line
    # ax.plot([np.pi/2, np.pi/2], [min_val, max_val], color='cyan', linewidth=5, alpha=.5, linestyle='-')
    # ax.plot([np.pi*3/2, np.pi*3/2], [min_val, max_val], color='cyan', linewidth=5, alpha=.5, linestyle='-')

    ## plot pref ori (noad) as a line
    ax.plot([ori_pref_noad, ori_pref_noad], [min_val, max_val], color='blue', linewidth=5, alpha=.5, label='pref ori no adapter')
    ori_pref_noad_opp = ori_pref_noad + np.pi # find opposite orientation of pref ori
    if ori_pref_noad_opp > 2*np.pi:
        ori_pref_noad_opp = ori_pref_noad_opp - 2*np.pi
    ax.plot([ori_pref_noad_opp, ori_pref_noad_opp], [min_val, max_val], color='blue', linewidth=5, alpha=.5, linestyle='-')

    ## plot pref ori (ad 250) as a line
    ax.plot([ori_pref_ad_250, ori_pref_ad_250], [min_val, max_val], color='orange', linewidth=5, alpha=.5, label='pref ori isi=250')
    ori_pref_ad_250_opp = ori_pref_ad_250 + np.pi
    if ori_pref_ad_250_opp > 2*np.pi:
        ori_pref_ad_250_opp = ori_pref_ad_250_opp - 2*np.pi
    ax.plot([ori_pref_ad_250_opp, ori_pref_ad_250_opp], [min_val, max_val], color='orange', linewidth=5, alpha=.5, linestyle='-')

    ## add text at top left
    bias_color = 'blue'
    minus_flag = 'pos'
    if tuning_bias < 0: # attractive bias
        bias_color = 'red'
        minus_flag = 'neg'
    plt.text(0.05, 0.9, f'tuning_bias: {np.round(tuning_bias, 2)}', transform=plt.gcf().transFigure, color=bias_color, fontsize=16) # in cartesian coordinates
    gOSI_color = 'green'
    if gOSI_noad < 0.5 or gOSI_250 < 0.5:
        gOSI_color = 'orange'
    if gOSI_noad < 0.5 and gOSI_250 < 0.5:
        gOSI_color = 'red'
    plt.text(0.05, 0.85, f'gOSI_noad: {np.round(gOSI_noad, 2)}', transform=plt.gcf().transFigure, color=gOSI_color, fontsize=16) # in cartesian coordinates
    plt.text(0.05, 0.8, f'gOSI_250: {np.round(gOSI_250, 2)}', transform=plt.gcf().transFigure, color=gOSI_color, fontsize=16) # in cartesian coordinates

    rticks = np.arange(0, max_val, step=0.05)
    ax.set_rticks(rticks)  # fewer radial ticks

    ax.set_xticks(np.linspace(0, 2*np.pi, 17))
    ax.set_xticklabels(np.arange(0, 360+22.5, 22.5))
    xticklabels = [label.get_text() for label in ax.get_xticklabels()]
    xticklabels[-1] = '' # set final xticklabel invisible
    ax.set_xticklabels(xticklabels)
    ax.set_ylim(min_val, max_val)

    area = df_tuning[df_tuning.cell_id == icell].area.values[0]
    ax.set_title(f'cell_id {icell} in {area}')
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.legend(bbox_to_anchor=(1.1, 1.1), frameon=False); # set legend position out of the way

    # break

    # save figure
    dir_result = r'C:\Users\lan\Documents\repos\inter\results\tuning bias single cell'.replace('\\', '/')
    fig.savefig(os.path.join(dir_result, f'gOSI_{np.round(gOSI_noad, 2)}_{icell}_tuning_curve.pdf'), bbox_inches='tight')
    plt.close(fig)
    clear_output(wait=True)


# ## filter cell gOSI
# for well fit cells, calculate global orientation selectivity index (gOSI)  
# formula taken from [Causal importance of orientation selectivity for generalization in image recognition](https://openreview.net/pdf?id=Bkx_Dj09tQ)

# In[156]:


for icell in tqdm(df_tuning.cell_id.unique()):
    df_cell = df_tuning[df_tuning.cell_id == icell]

    tuning_noad = df_cell.tuning_noad.values[0] - min(df_cell.tuning_noad.values[0]) # ensure all values are non negative
    tuning_250 = df_cell.tuning_250.values[0] - min(df_cell.tuning_250.values[0])

    theta_arr = np.linspace(0, 180-22.5, 8) # according to formula: unit deg, not rad
    sin_arr = np.sin(2 * theta_arr)
    cos_arr = np.cos(2 * theta_arr)

    gOSI_noad = np.sqrt((np.sum(tuning_noad * sin_arr))**2 + (np.sum(tuning_noad * cos_arr))**2) / np.sum(tuning_noad)
    gOSI_250 = np.sqrt((np.sum(tuning_250 * sin_arr))**2 + (np.sum(tuning_250 * cos_arr))**2) / np.sum(tuning_250)

    df_tuning.loc[df_tuning.cell_id == icell, 'gOSI_noad'] = gOSI_noad
    df_tuning.loc[df_tuning.cell_id == icell, 'gOSI_250'] = gOSI_250

df_tuning = df_tuning.sort_values(by=['gOSI_noad', 'gOSI_250'], ascending=False)
df_tuning.groupby('cell_id', sort=False)['area', 'gOSI_noad', 'gOSI_250'].first()#.head(20)


# In[26]:


# fig, axes = plt.subplots(1, 2, figsize=(10, 5))
# df_tuning.sort_values(by='area', ascending=False).groupby(['area', 'cell_id'], sort=False).gOSI_noad.hist(bins=20, alpha=.5, ax=axes[0]);
# df_tuning.sort_values(by='area', ascending=False).groupby(['area', 'cell_id'], sort=False).gOSI_250.hist(bins=20, alpha=.5, ax=axes[1]);

# axes[0].legend(['V1', 'LM', 'LI',], frameon=False) # sorted area descending
# plt.tight_layout();
# ## indeed, gOSI is between 0 and 1. higher gOSI means more orientation selective


# In[27]:


df_tuning.groupby('area').gOSI_noad.describe() # df tuning is well-fit cells only. even so, LI gOSI is much worse than V1 and LM


# In[28]:


df_tidy.groupby('area').filter_cell_well_fit.describe() # well fit cell filter discards most cells in LI already


# In[29]:


tmp1 = df_tuning.sort_values(by=['gOSI_noad', 'gOSI_250'], ascending=False).gOSI_noad.values
tmp2 = df_tuning.sort_values(by=['gOSI_noad', 'gOSI_250'], ascending=False).gOSI_250.values
tmp3 = df_tuning.sort_values(by=['gOSI_noad', 'gOSI_250'], ascending=False).tuning_bias.values

plt.scatter(tmp1, tmp2, alpha=.1, s=1);
plt.plot([0, 1], [0, 1], color='k', linestyle='-', linewidth=1); # draw diagonal line

plt.xlim([0, 1]);
plt.ylim([0, 1]);
plt.xticks([0, .5, 1]);
plt.yticks([0, .5, 1]);
plt.xlabel('gOSI_noad');
plt.ylabel('gOSI_250');
plt.gca().set_aspect('equal', adjustable='box'); # set axis square

## gOSI noad is usually higher. gOSI noad vs 250 is correlated


# In[30]:


## regression: gOSI vs tuning_bias
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score

tmp1 = df_tuning.sort_values(by=['gOSI_noad', 'gOSI_250'], ascending=False).gOSI_noad.values
tmp2 = df_tuning.sort_values(by=['gOSI_noad', 'gOSI_250'], ascending=False).gOSI_250.values
tmp3 = df_tuning.sort_values(by=['gOSI_noad', 'gOSI_250'], ascending=False).tuning_bias.values
tmp1 = tmp1[~np.isnan(tmp3)]
tmp2 = tmp2[~np.isnan(tmp3)]
tmp3 = tmp3[~np.isnan(tmp3)]
tmp1 = tmp1.reshape(-1, 1)
tmp2 = tmp2.reshape(-1, 1)
tmp3 = tmp3.reshape(-1, 1)

reg1 = LinearRegression().fit(tmp1, tmp3)
reg2 = LinearRegression().fit(tmp2, tmp3)
print('r2 score of gOSI_noad vs tuning_bias: ', r2_score(tmp3, reg1.predict(tmp1)))
print('r2 score of gOSI_250 vs tuning_bias: ', r2_score(tmp3, reg2.predict(tmp2)))

plt.scatter(tmp1, tmp3, alpha=.1, s=1);
plt.plot(tmp1, reg1.predict(tmp1), color='k', linestyle='-', linewidth=1);
print('slope: ', reg1.coef_[0][0])
print('intercept: ', reg1.intercept_[0])

print('positive correlation between gOSI and tuning_bias')


# ### tuning bias plot - high gOSI
# didnt seem to help

# In[265]:


gOSI_thres = 0.5 # In visual neuroscience, neurons with gOSI > 0.33 are often considered to be orientation-selective (Piscopo et al., 2013; Kondo & Ohki, 2015).
df_tuning_gOSI = df_tuning[(df_tuning.gOSI_noad > gOSI_thres) & (df_tuning.gOSI_250 > gOSI_thres)]
print(df_tuning_gOSI.groupby('area').cell_id.nunique())

df_tuning_gOSI = df_tuning_gOSI.groupby('cell_id')[['area', 'pref_unadapted_distance_bin', 'tuning_bias', 'tuning_bias_control']].first().reset_index()
df_tuning_gOSI


# In[266]:


fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
g1 = sns.pointplot(x='pref_unadapted_distance_bin', y='tuning_bias', 
                   data=df_tuning_gOSI[df_tuning_gOSI.area != 'LI'], hue='area',
                   errorbar='sd', # ('ci', 68), 
                   errwidthfloat=1, capsize=.1,
                   ax=axes[0], dodge=True,
                   )
g2 = sns.pointplot(x='pref_unadapted_distance_bin', y='tuning_bias_control', 
                   data=df_tuning_gOSI[df_tuning_gOSI.area != 'LI'], hue='area',
                  errorbar=('ci', 68), errwidthfloat=1, capsize=.1,
                  ax=axes[1], dodge=True, legend=False,
                  )
axes[0].axhline(0, color='gray', linestyle='-', alpha=.3);
axes[1].axhline(0, color='gray', linestyle='-', alpha=.3);

g2.set(ylabel=None) # remove ylabel
plt.setp(axes[0].collections, alpha=.5);
plt.setp(axes[1].collections, alpha=.5);
plt.setp(axes[0].lines, alpha=.7);
plt.setp(axes[1].lines, alpha=.7);
fig.tight_layout();

# plt.savefig(os.path.join(dir_fig, 'tuning_bias_by_area.pdf'))


# In[267]:


fig, axes = plt.subplots(1, 2, figsize=(12, 5), sharey=True)
g1 = sns.pointplot(x='pref_unadapted_distance_bin', y='tuning_bias', data=df_tuning_gOSI, hue='area',
                   errorbar='sd', # ('ci', 68), 
                   errwidthfloat=1, capsize=.1,
                   ax=axes[0], dodge=True,
                   )
g2 = sns.pointplot(x='pref_unadapted_distance_bin', y='tuning_bias_control', data=df_tuning_gOSI, hue='area',
                  errorbar=('ci', 68), errwidthfloat=1, capsize=.1,
                  ax=axes[1], dodge=True, legend=False,
                  )
axes[0].axhline(0, color='gray', linestyle='-', alpha=.3);
axes[1].axhline(0, color='gray', linestyle='-', alpha=.3);

g2.set(ylabel=None) # remove ylabel
plt.setp(axes[0].collections, alpha=.5);
plt.setp(axes[1].collections, alpha=.5);
plt.setp(axes[0].lines, alpha=.7);
plt.setp(axes[1].lines, alpha=.7);
fig.tight_layout();

# plt.savefig(os.path.join(dir_fig, 'tuning_bias_by_area.pdf'))


# In[268]:


df_tuning_gOSI['tuning_bias_diff'] = df_tuning_gOSI['tuning_bias'] - df_tuning_gOSI['tuning_bias_control']

fig, axes = plt.subplots(1, 1, figsize=(9, 5), sharey=True)
g1 = sns.pointplot(x='pref_unadapted_distance_bin', y='tuning_bias_diff', data=df_tuning_gOSI, hue='area',
                   errorbar=('ci', 68), errwidthfloat=1, capsize=.1,
                   ax=axes, dodge=True,
                   )

# annotation above each dot, ncell
ncell_bin = df_tuning_gOSI.groupby('pref_unadapted_distance_bin').cell_id.nunique().sort_index().values
ylim = axes.get_ylim()
axes.set_ylim(ylim[0], ylim[1] + 1.5)
for i in range(len(ncell_bin)):
    axes.annotate(str(ncell_bin[i]), (i, ylim[1] + 0.1), ha='center', va='center', size=12)

axes.legend(frameon=False, bbox_to_anchor=(1.5, 0.5))
g1.set(xlabel='pref ori distance from 0 deg adapter')
g1.set(ylabel='tuning bias diff (real - control)')
axes.axhline(0, color='gray', linestyle='-', alpha=.3);
plt.setp(axes.collections, alpha=.5);
plt.setp(axes.lines, alpha=.7);
fig.tight_layout();

# plt.savefig(os.path.join(dir_fig, 'tuning_bias_control_diff_by_area.pdf'))


# ## filter cell ori-mod
# dont fit tuning curve, dont filter well fit cells.  
# find orientation-modulated cells by anova across all orientations responses

# In[31]:


## construct tuning_vec column

cell_property = (df_tidy[['cell_id', 'isi', 'area', 'filter_cell_vis']] # need cell info: area, vis driven
                 .groupby(['cell_id', 'isi']) # prepare to match with df_ori_mod
                 .first() # only take first value. all values should be the same for each cell and isi combination
                 .reset_index())

df_ori_mod = (df_tidy[df_tidy.resp_id == 'R2'] # only R2 has diff ori
            [['dfof', 'cell_id', 'resp_id', 'isi', 'stim2_id']]
              .groupby(['cell_id', 'isi', 'stim2_id']).agg({'dfof': 'mean'}) # aggregate resp by cell, isi, ori
              .groupby(['cell_id', 'isi']).agg({'dfof':lambda x: list(x)}) # each row is a list of aggregated resp across ori
              .reset_index()
              .rename(columns={'dfof': 'tuning_vec'})
              .merge(cell_property, on=['cell_id', 'isi'], how='left') # merge with cell info
            )
df_ori_mod


# In[32]:


## construct ori_mod column

from scipy.stats import kruskal
p_threshold = 0.05

df_kruskal = (df_tidy[df_tidy.resp_id == 'R2'] # only R2 has diff ori
            [['dfof', 'cell_id', 'isi', 'stim2_id']]
            .groupby(['cell_id', 'isi', 'stim2_id']) # for each cell, each isi condition, calc responses to each ori
            .agg({'dfof':lambda x: list(x)}) # each row, dfof col contains a list of ori responses across trials
            .reset_index()
            .groupby(['cell_id', 'isi']) # for each cell, each isi condition, calc ori modulation
            .apply(lambda x: kruskal(*x.dfof.values).pvalue < p_threshold) # kruskal, where each ori is a group
            .reset_index()
            )
df_kruskal = df_kruskal.rename(columns={0: 'ori_mod'}) # bool col for ori modulation

df_ori_mod = df_ori_mod.merge(df_kruskal, on=['cell_id', 'isi'], how='left')
df_ori_mod


# In[33]:


## construct max_ori and max_ori_distance column

df_ori_mod['max_ori'] = df_ori_mod.tuning_vec.apply(lambda x: np.argmax(x)) # NOTE: even though we calculated max_ori for isi 250 and 750, only isi 6000 was actually used below (inherited by df_lineplot)
df_ori_mod['max_ori_dist'] = df_ori_mod.max_ori.apply(lambda x: 22.5*(8-x) if x > 4 else 22.5*x) # distance from 0 deg, hard coded for 8 ori. TODO: modify if needed

# df_ori_mod.groupby('max_ori').max_ori_distance.value_counts()
# df_ori_mod.groupby('max_ori_distance').max_ori.value_counts()
df_ori_mod.groupby('max_ori').max_ori_dist.unique()


# In[34]:


## prepare for lineplot. don't need tuning vec, but need (stim2_id, dfof) observations
df_lineplot = (df_tidy[df_tidy.resp_id == 'R2'] # only R2 has diff ori
                      [['dfof', 'cell_id', 
                        'area', 'filter_cell_vis',
                        'resp_id', 'isi', 'stim2_id']]
                        .reset_index(drop=True)
                        )
## inherit ori_mod, max_ori, max_ori_distance from df_ori_mod
df_lineplot = df_lineplot.merge(
      df_ori_mod[df_ori_mod.isi == 6000]
      [['cell_id', 'ori_mod', 'max_ori', 'max_ori_dist']], 
      on=['cell_id'], how='left') # only use isi 6000 (no adapter condition) to determine ori_mod, max_ori and max_ori_dist for each cell

## inherit tuning_vec from df_ori_mod, for each isi
df_lineplot = df_lineplot.merge(
      df_ori_mod
      [['cell_id', 'tuning_vec', 'isi']],
      on=['cell_id', 'isi'], how='left') # take tuning_vec from df_ori_mod, for each cell and isi

df_lineplot


# In[35]:


# for each cell_id, divide tuning_vec of isi 250 by tuning_vec of isi 6000
tmp = df_lineplot.groupby(['cell_id', 'isi']).tuning_vec.first().to_frame()
tmp

## split tuning_vec into 8 columns
tmp = tmp.tuning_vec.apply(pd.Series)
tmp

tmp_noad = (tmp.groupby(level=['cell_id']).transform('last'))
tmp_noad

result = tmp - tmp_noad
result

## merge 8 columns back to tuning_vec
result['tuning_diff']= result.values.tolist()
result


# In[36]:


plt.hist(tmp_noad.values.flatten(), bins=100);
plt.xlim([-0.1, 0.25])

np.min(tmp_noad.values.flatten()), np.max(tmp_noad.values.flatten()), \
    sum(tmp_noad.values.flatten() < 0) / len(tmp_noad.values.flatten())


# ## tuning curve bias
# no fitting, just align max response ori

# In[37]:


## merge tuning_diff back to df_lineplot
df_lineplot = df_lineplot.merge(
        result.loc[:, 'tuning_diff'].reset_index(), # NOTE tuning_diff is (tuning_ad250 - tuning_noad)
        on=['cell_id', 'isi'], how='left')
df_lineplot


# In[49]:


df_lineplot['date'] = df_lineplot.cell_id.apply(lambda x: x.split('_')[0])
df_lineplot['sess'] = df_lineplot.cell_id.apply(lambda x: x.split('_')[1])
df_lineplot[df_lineplot.area == 'V1'].date.unique()


# In[54]:


df_filter = df_lineplot[(df_lineplot.filter_cell_vis == True)
                       & (df_lineplot.ori_mod == True)
                    #    & (df_lineplot.isi != 750)
                       & (df_lineplot.area == 'LM')
                     ##  & (df_lineplot.date == '201015') # temporary filter to date, to match matlab san check for TC
                      #  & (df_lineplot.date == '210120') # temporary filter to date, to match matlab san check for TC
                       ]
# print(df_filter.groupby('max_ori').cell_id.nunique())
# 
df_tuning_diff = df_filter.groupby(['cell_id', 'isi']).tuning_diff.first().to_frame()
df_tuning_diff = df_tuning_diff.tuning_diff.apply(pd.Series)

# rename columns by adding 'stim' to each column
df_tuning_diff = df_tuning_diff.rename(columns={i: 'stim' + str(i) for i in range(8)})
df_tuning_diff = df_tuning_diff.reset_index()
df_tuning_diff

# pivot table so stim is row
df_tuning_diff = df_tuning_diff.melt(id_vars=['cell_id', 'isi'], var_name='stim', value_name='tuning_diff')
df_tuning_diff

# for stim col, take only the last char
df_tuning_diff['stim'] = df_tuning_diff.stim.apply(lambda x: x[-1])
df_tuning_diff.stim = df_tuning_diff.stim.astype(int)
df_tuning_diff

## set sns color palette
sns.lineplot(data=df_tuning_diff, x='stim', y='tuning_diff', hue='isi', ci=68, #alpha=0.9,
            # estimator=np.median, 
            palette=['r', 'orange', 'b'],
            )
plt.xticks([0, 4], ['0', '90'])
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# In[55]:


df_filter = df_lineplot[(df_lineplot.area == 'V1') # check V1 only
                       & (df_lineplot.ori_mod == True)
                       & (df_lineplot.filter_cell_vis == True)
                       ]

sns.lineplot(data=df_filter, x='stim2_id', y='dfof', hue='isi', ci=68, alpha=0.6, 
            estimator=np.median, 
            palette=['r', 'orange', 'b'],)
plt.xticks([0, 4], ['0', '90'])
plt.xlabel('stim2 orientation (deg)')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# In[56]:


sns.lineplot(data=df_filter, x='max_ori_dist', y='dfof', hue='isi', ci=68, alpha=0.6, 
            estimator=np.mean, 
            palette=['r', 'orange', 'b'],)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# In[47]:


df_filter = df_lineplot[(df_lineplot.filter_cell_vis == True)
                       & (df_lineplot.ori_mod == True)
                       & (df_lineplot.isi != 750)
                       & (df_lineplot.area == 'V1')
                     ##  & (df_lineplot.date == '201015') # temporary filter to date, to match matlab san check for TC
                      #  & (df_lineplot.date == '210120') # temporary filter to date, to match matlab san check for TC
                       ]

sns.set(font_scale=0.8, context='talk', style='whitegrid')
sns.lineplot(data=df_filter, x='stim2_id', y='dfof', hue='isi', ci=68, alpha=0.6, 
            # estimator=np.median
            )
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);

# sns.set(font_scale=0.8, context='talk', style='whitegrid')
# g = sns.FacetGrid(df_filter, row="date", hue="isi")
# g.map_dataframe(sns.lineplot, x="stim2_id", y="dfof", alpha=0.6)
# g.figure.subplots_adjust(wspace=0.1, hspace=0.3)
# g.add_legend();

# # set xticks to 0, 4
# for ax in g.axes.flat:
#     ax.set_xticks([0, 4])
#     ax.set_xticklabels(['0', '90'])


# In[394]:


df_filter = df_lineplot[(df_lineplot.filter_cell_vis == True)
                       & (df_lineplot.ori_mod == True)
                       & (df_lineplot.isi != 750)
                       & (df_lineplot.area == 'LM')
                     ##  & (df_lineplot.date == '201015') # temporary filter to date, to match matlab san check for TC
                      #  & (df_lineplot.date == '210120') # temporary filter to date, to match matlab san check for TC
                       ]

sns.set(font_scale=0.8, context='talk', style='whitegrid')
sns.lineplot(data=df_filter, x='stim2_id', y='dfof', hue='isi', ci=68, alpha=0.6, estimator=np.median)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# In[395]:


df_filter = df_lineplot[(df_lineplot.filter_cell_vis == True)
                       & (df_lineplot.ori_mod == True)
                       & (df_lineplot.isi != 750)
                       & (df_lineplot.area == 'LI')
                     ##  & (df_lineplot.date == '201015') # temporary filter to date, to match matlab san check for TC
                      #  & (df_lineplot.date == '210120') # temporary filter to date, to match matlab san check for TC
                       ]

sns.set(font_scale=0.8, context='talk', style='whitegrid')
sns.lineplot(data=df_filter, x='stim2_id', y='dfof', hue='isi', ci=68, alpha=0.6, estimator=np.median)
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# ## Jin 2019 Fig 2

# In[66]:


df_lineplot.area.unique()


# In[76]:


## construct df for tuning_diff lineplot:
## for each cell_id, subtract tuning_vec of isi 250 by tuning_vec of isi 6000
tmp = df_lineplot.groupby(['cell_id', 'isi']).tuning_vec.first().to_frame()
tmp.head(6)

tmp_noad_shift = tmp.copy()
tmp_noad_shift['tuning_vec'] = tmp_noad_shift.tuning_vec.apply(lambda x: [min(x) for i in range(len(x))])
tmp_noad_shift.head(6)

tmp_noad_shift = (tmp_noad_shift.groupby(level=['cell_id']).transform('last')) # last isi is 6000
epsilon = 0.001
tmp_noad_shift.tuning_vec = tmp_noad_shift.tuning_vec.apply(lambda x: [-min(0, i) + epsilon for i in x])
tmp_noad_shift.head(6)

## split tuning_vec into 8 columns
tmp = tmp.tuning_vec.apply(pd.Series)
tmp.head(6)
# plt.hist(tmp.values.flatten(), bins=100);
# min(tmp.values.flatten()), max(tmp.values.flatten())

tmp_noad_shift = tmp_noad_shift.tuning_vec.apply(pd.Series)
tmp_noad_shift.head(6)

tmp = tmp + tmp_noad_shift # shift tuning_vec of all isi by min(tuning_vec) of isi 6000
tmp.head(6)

tmp_noad = (tmp.groupby(level=['cell_id']).transform('last'))
tmp_noad.head(6)

result = tmp / tmp_noad
result
# plt.hist(result.values.flatten(), bins=100);
# min(result.values.flatten()), max(result.values.flatten())

thresh = 5
result = result.applymap(lambda x: x if abs(x) < thresh else np.nan) # threshold too large values - change to nan
result

## merge 8 columns back to tuning_vec
result['tuning_diff']= result.values.tolist()
result

## merge tuning_diff back to df_lineplot
df_lineplot = df_lineplot.loc[:, df_lineplot.columns != 'tuning_diff'].merge( # prevent duplicate column
        result.loc[:, 'tuning_diff'].reset_index(), # NOTE tuning_diff is (tuning_ad250 - tuning_noad)
        on=['cell_id', 'isi'], how='left')
df_lineplot


# In[94]:


df_filter = df_lineplot[(df_lineplot.area == 'V1')
                       & (df_lineplot.ori_mod == True)
                       & (df_lineplot.filter_cell_vis == True)
                       ]
# df_filter

df_tuning_diff = df_filter.groupby(['cell_id', 'isi']).tuning_diff.first().to_frame()
df_tuning_diff = df_tuning_diff.tuning_diff.apply(pd.Series)

# rename columns by adding 'stim' to each column
df_tuning_diff = df_tuning_diff.rename(columns={i: 'stim' + str(i) for i in range(8)})
df_tuning_diff = df_tuning_diff.reset_index()
df_tuning_diff

# pivot table so stim is row
df_tuning_diff = df_tuning_diff.melt(id_vars=['cell_id', 'isi'], var_name='stim', value_name='tuning_diff')
df_tuning_diff

# for stim col, take only the last char
df_tuning_diff['stim'] = df_tuning_diff.stim.apply(lambda x: x[-1])
df_tuning_diff.stim = df_tuning_diff.stim.astype(int)
df_tuning_diff['ori_dist'] = df_tuning_diff.stim.apply(lambda x: 22.5*(8-x) if x > 4 else 22.5*x) # stim2 distance from 0 deg, hard coded for 8 ori. TODO: modify if needed. here, we dont align by max ori of cells (to group cells), instead we plot all cells
df_tuning_diff

## inherit max_ori from df_lineplot
df_tuning_diff = df_tuning_diff.merge(df_filter[['cell_id', 'isi', 'max_ori', 'max_ori_dist']], on=['cell_id', 'isi'], how='inner')
df_tuning_diff = df_tuning_diff.drop_duplicates() # TODO: why are there duplicates?
df_tuning_diff


# In[95]:


df_tuning_diff.cell_id.nunique() # how many cells


# In[96]:


## map with dict
ori_dist_bin = {0:0, 
                22.5:45,
                45:45,
                67.5:90,
                90:90,}

df_tuning_diff['ori_dist_bin'] = df_tuning_diff.ori_dist.map(ori_dist_bin)
df_tuning_diff['max_ori_dist_bin'] = df_tuning_diff.max_ori_dist.map(ori_dist_bin)

df_tuning_diff


# In[97]:


tmp = df_tuning_diff[df_tuning_diff.max_ori == df_tuning_diff.stim] # only take rows with max ori = stim2 ori
tmp

sns.lineplot(data=tmp, x='max_ori_dist_bin', y='tuning_diff', hue='isi', ci=68, #alpha=0.9,
            # estimator=np.median, 
            palette=['r', 'orange', 'b'],
            )

# set xticks 0, 45, 90
plt.xticks([0, 45, 90], [0, 45, 90])

plt.xlabel('Pref distance from adapter (deg)')
plt.ylabel('Tuning after adaptation / before')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);

# dir_fig = r'C:\Users\ll357\Documents\inter\results\tuning curve bias san check/'.replace('\\', '/')
# plt.savefig(dir_fig + 'Fig2E Jin2019.pdf', bbox_inches='tight')


# In[98]:


sns.lineplot(data=tmp, x='max_ori_dist', y='tuning_diff', hue='isi', ci=68, #alpha=0.9,
            # estimator=np.median, 
            palette=['r', 'orange', 'b'],
            )
plt.xlabel('Pref distance from adapter (deg)')
plt.ylabel('Tuning after adaptation / before')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);

# dir_fig = r'C:\Users\ll357\Documents\inter\results\tuning curve bias san check/'.replace('\\', '/')
# plt.savefig(dir_fig + 'tuning_diff_lineplot_binned.pdf', bbox_inches='tight')


# In[89]:


sns.lineplot(data=df_tuning_diff, x='ori_dist_bin', y='tuning_diff', hue='isi', ci=68, #alpha=0.9,
            # estimator=np.median, 
            palette=['r', 'orange', 'b'],
            )
plt.xlabel('Stim2 distance from adapter (deg)')
plt.ylabel('Tuning after adaptation / before')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);

# dir_fig = r'C:\Users\ll357\Documents\inter\results\tuning curve bias san check/'.replace('\\', '/')
# plt.savefig(dir_fig + 'tuning vs stim2_binned-adapter dist.pdf', bbox_inches='tight')


# In[90]:


## set sns color palette
sns.lineplot(data=df_tuning_diff, x='ori_dist', y='tuning_diff', hue='isi', ci=68, #alpha=0.9,
            # estimator=np.median, 
            palette=['r', 'orange', 'b'],
            )
plt.xlabel('Stim2 distance from adapter (deg)')
plt.ylabel('Tuning after adaptation / before')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);

# dir_fig = r'C:\Users\ll357\Documents\inter\results\tuning curve bias san check/'.replace('\\', '/')
# # plt.savefig(dir_fig + 'tuning vs stim2-adapter dist.pdf', bbox_inches='tight')
# plt.savefig(dir_fig + 'Fig2D Jin2019.pdf', bbox_inches='tight')


# ## facet grid

# In[378]:


df_filter = df_lineplot[(df_lineplot.filter_cell_vis == True)
                      #  & (df_lineplot.ori_mod == True)
                       & (df_lineplot.isi != 750)
                       & (df_lineplot.area == 'V1')
                     ##  & (df_lineplot.date == '201015') # temporary filter to date, to match matlab san check for TC
                      #  & (df_lineplot.date == '210120') # temporary filter to date, to match matlab san check for TC
                       ]

sns.set(font_scale=0.8, context='talk', style='whitegrid')
# g = sns.FacetGrid(df_filter, col="area", row="max_ori_dist", hue="isi")
g = sns.FacetGrid(df_filter, col="max_ori_dist", row="date", hue="isi")
# g = sns.FacetGrid(df_filter, col="date", row="max_ori_dist", hue="isi") # which date has the most outrageous adp at orthogonal ori relative to adapter ori
g.map_dataframe(sns.lineplot, x="stim2_id", y="dfof", alpha=0.6)
g.figure.subplots_adjust(wspace=0.45, hspace=0.3)
g.add_legend();

# set xticks to 0, 4
for ax in g.axes.flat:
    ax.set_xticks([0, 4])
    ax.set_xticklabels(['0', '90'])


# In[379]:


sns.set(font_scale=0.8, context='talk', style='whitegrid')
# g = sns.FacetGrid(df_filter, col="area", row="max_ori", hue="isi")
g = sns.FacetGrid(df_filter, col="date", row="max_ori", hue="isi")
g.map_dataframe(sns.lineplot, x="stim2_id", y="dfof", alpha=0.6)
g.figure.subplots_adjust(wspace=0.4, hspace=0.3)
g.add_legend();

for ax in g.axes.flat:
    ax.set_xticks([0, 4])
    ax.set_xticklabels(['0', '90'])

dir_result = r'C:\Users\ll357\Documents\inter\results\tuning curve bias san check'.replace('\\', os.sep)
g.savefig(os.path.join(dir_result, 'tuning_before_vs_after_adp.pdf'), bbox_inches='tight')


# # san check
# what can generate possibly fake adp when adapter vs target are orthogonal?  
# check if bin=90 tuning curve is real: plot timecourse for stim2=90, noad vs 250

# In[335]:


df_test = df_lineplot[(df_lineplot.date == '210120')] # use this date and sess for sanity check
df_test.cell_id = df_test.cell_id.apply(lambda x: x.split('_')[-1]) # temporarily remove date from cell_id
df_test


# In[369]:


df_filter = df_test[(df_test.filter_cell_vis == True)
                    # & (df_test.ori_mod == True)
                    & (df_test.isi != 750)
                    ]
df_filter.cell_id.nunique(), df_test.cell_id.nunique()

cell_id_bool = pd.DataFrame(df_test.cell_id.astype(int).unique()).isin(df_filter.cell_id.astype(int).unique())
cell_id_bool.head(20) # whether this cell passes the filter (vis and ori_mod)
cell_id_bool = cell_id_bool.values.astype(int).flatten()
cell_id_bool

import scipy.io as sio
dir_mat = r'C:\Users\ll357\Documents\inter\results\tuning curve bias san check/'.replace('\\', '/')
# sio.savemat(dir_mat+'vis_orimod_cell_bool.mat', {'vis_orimod_cell_bool': cell_id_bool})
sio.savemat(dir_mat+'vis_cell_bool.mat', {'vis_cell_bool': cell_id_bool}) # if only filter by vis, not ori_mod


# # decorrelation
# cosine similarity

# In[1052]:


## suppress SettingWithCopyWarning
pd.options.mode.chained_assignment = None  # default='warn'

df_decorr = df_tidy[['dfof', 
                    'cell_id', 'resp_id', 'isi', 'stim2_id', 
                    'filter_cell_vis', 'filter_cell_stim', 'filter_cell_well_fit', 
                    'area', 'mouse', 'date', 'sess'
                    ]]

df_decorr['stim2_dist'] = df_decorr.stim2_id.apply(lambda x: 22.5*(8-x) if x > 4 else 22.5*x) 
                          # stim2 distance from 0 deg, hard coded for 8 ori. TODO: modify if needed
df_decorr['date_sess'] = df_decorr.date + '_' + df_decorr.sess

current_area = 'V1'
df_decorr = df_decorr[(df_decorr.isi != 750)
                      & (df_decorr.filter_cell_vis == True)
                      # & (df_decorr.filter_cell_stim == True)
                      # & (df_decorr.filter_cell_well_fit == True)
                      & (df_decorr.resp_id == 'R2')
                      & (df_decorr.area == current_area)
                    ]

df_decorr


# ## decorr from self
# before vs after adp

# In[467]:


cos_sim_sets = []
for iset in tqdm(df_decorr.date_sess.unique()):
    # print(iset)
    
    df_pop_vec = (df_decorr[df_decorr.date_sess == iset]
                  .groupby(['isi', 'stim2_dist', 'cell_id'])
                  .dfof.mean().reset_index()
                  .pivot_table(index=['isi', 'stim2_dist'], columns='cell_id', values='dfof')
                )
    nori = df_decorr.stim2_dist.nunique()
    
    cos_sim_iset = []
    for iori in range(nori):
        pop_vec_ref = df_pop_vec.iloc[nori + iori, :].values # ref is the same ori, but before adp (isi=6000)
        pop_vec_ori = df_pop_vec.iloc[iori, :].values # pop vec after adp (isi=250)
        cos_sim = dot(pop_vec_ref, pop_vec_ori) / (norm(pop_vec_ref) * norm(pop_vec_ori))
        cos_sim_iset.append(cos_sim)
    cos_sim_sets.append(cos_sim_iset)
    
len(cos_sim_sets), len(cos_sim_sets[0]) # n of recording, n of self-decorrelation pairs (nori)


# In[468]:


cos_sim_sets = np.array(cos_sim_sets)
cos_sim_sets_avg = np.nanmean(cos_sim_sets, axis=0) # across recordings
cos_sim_sets_sem = np.nanstd(cos_sim_sets, axis=0) / np.sqrt(cos_sim_sets.shape[0])

fig, ax = plt.subplots(1, 1, figsize=(6, 4), sharex=True)
ax.errorbar(x=np.arange(len(cos_sim_sets_avg)),
               y=cos_sim_sets_avg,
               yerr=cos_sim_sets_sem,
               color='blue', # label='isi 250', 
               alpha=0.7)

# ax.legend(frameon=False, loc='lower center')
ax.set_xlabel('stim2 distance from adapter (deg)');
ax.set_ylabel('cosine similarity')
ax.set_xticks(range(nori), ['0', '22.5', '45', '67.5', '90'])

fig.tight_layout();
fig_dir = r'C:\Users\lan\Documents\repos\inter\results\decorrelation vs adp\ref=self before adp'.replace('\\', os.sep)
fig_name = 'decorr_from_self_after_adp_' + current_area + '.pdf'
plt.savefig(os.path.join(fig_dir, fig_name), bbox_inches='tight')


# ## decorr from ref repre

# In[1044]:


(df_tidy[(df_tidy.filter_cell_vis == True)
        #  & (df_tidy.filter_cell_well_fit == True)
         ]
        .groupby('area')
        .cell_id.nunique())


# In[ ]:


''' 
if resample cells without affecting trials: 
1. groupby cell_id, then sample n cells. set cell_id as index, use df.reindex(cell_id) to get the resampled df
2. set cell_id as index in df_decorr, df_decorr[df_decorr.cell_id.isin(subsample_cell)] -> df.reindex
3. df_pop_vec. ... .groupby, reset_index, sample n cells, pivot_table
'''


# In[1053]:


## bootstrap resample cells to get error bar

nboot = 100
# subsample_cell_ratio = 0.7
# ncell_subsample = (df_tidy[(df_tidy.filter_cell_vis == True)]
#                    .groupby('area')
#                    .cell_id.nunique()
#                    .min()) # min ncells across areas - determine subsample size
cos_sim_boots = []

for iboot in tqdm(range(nboot)):
    
    np.random.seed(iboot)
    subsample_cell = np.random.choice(df_decorr.cell_id.unique(), 
                                    #   size=int(subsample_cell_ratio*df_decorr.cell_id.nunique()), 
                                    #   size=ncell_subsample,
                                      size=df_decorr.cell_id.nunique(), 
                                      replace=True)
    
    # df_pop_vec = df_decorr.sample(frac=1, 
    #                               replace=True, 
    #                               random_state=iboot) # resample cellxtrial (row-wise) with replacement, same size as original
    
    df_pop_vec = (df_decorr # [df_decorr.cell_id.isin(subsample_cell)]
                  .groupby(['isi', 'stim2_dist', 'cell_id'])
                  .dfof.median().reset_index() # cell-level median, under each isi-stim2_dist condition
                  .groupby(['isi', 'stim2_dist'])
                  .sample(frac=3, replace=True, random_state=iboot) # resample cells with replacement in each cond
                  .pivot_table(index=['isi', 'stim2_dist'], columns='cell_id', values='dfof') # df allows duplicate columns
                  .fillna(0) # fill nan with 0, bc some cells are not sampled in some conditions
                )
    nori = df_decorr.stim2_dist.nunique()
    nisi_now = df_decorr.isi.nunique() # discarded isi 750

    cos_sim_iboot = []
    for iisi in range(nisi_now):
        if iisi == 0:
            # iref = 0
            iref = 4
        elif iisi == 1:
            # iref = nori
            iref = nori+4 # -1 (0-based indexing) and +1 (go to next isi) cancel out
        pop_vec_ref = df_pop_vec.iloc[iref, :].values # reference pop vec: 0 deg target resp
        
        for iori in range(nori):
            irow = iisi * nori + iori
            pop_vec_ori = df_pop_vec.iloc[irow, :].values
            cos_sim = dot(pop_vec_ref, pop_vec_ori) / (norm(pop_vec_ref) * norm(pop_vec_ori))
            cos_sim_iboot.append(cos_sim)
    cos_sim_boots.append(cos_sim_iboot)
    
len(cos_sim_boots), len(cos_sim_boots[0]) # n of boot, n of decorrelation pairs (nisi * nori)


# In[1011]:


# tmp = (df_decorr # [df_decorr.cell_id.isin(subsample_cell)]
#         .groupby(['isi', 'stim2_dist', 'cell_id'])
#         .dfof.median().reset_index() # cell-level median, under each isi-stim2_dist condition
#         .groupby(['isi', 'stim2_dist'])
#         .sample(frac=4/122, replace=True, random_state=0) # resample cells with replacement, same size as original
#         .groupby(['isi', 'stim2_dist'])
#         # .pivot_table(index=['isi', 'stim2_dist'], columns='cell_id', values='dfof')
#         )
# tmp#.cell_id.nunique(), df_decorr.cell_id.nunique()

# for key, group in tmp:
#     # print(key)
#     # print(group)
#     break
# group.cell_id.nunique(), df_decorr.cell_id.nunique()
# group


# In[603]:


## separate recordings to get error bar

# cos_sim_sets = []
# for iset in tqdm(df_decorr.date_sess.unique()):
#     # print(iset)
    
#     df_pop_vec = (df_decorr[df_decorr.date_sess == iset]
#                   .groupby(['isi', 'stim2_dist', 'cell_id'])
#                   .dfof.mean().reset_index()
#                   .pivot_table(index=['isi', 'stim2_dist'], columns='cell_id', values='dfof')
#                 )
#     nori = df_decorr.stim2_dist.nunique()
#     nisi_now = df_decorr.isi.nunique() # discarded isi 750

#     cos_sim_iset = []
#     for iisi in range(nisi_now):
#         if iisi == 0:
#             iref = 0
#             # iref = 4
#         elif iisi == 1:
#             iref = nori
#             # iref = nori+4 # -1 (0-based indexing) and +1 (go to next isi) cancel out
#         pop_vec_ref = df_pop_vec.iloc[iref, :].values # reference pop vec: 0 deg target resp
        
#         for iori in range(nori):
#             irow = iisi * nori + iori
#             pop_vec_ori = df_pop_vec.iloc[irow, :].values
#             cos_sim = dot(pop_vec_ref, pop_vec_ori) / (norm(pop_vec_ref) * norm(pop_vec_ori))
#             cos_sim_iset.append(cos_sim)
#     cos_sim_sets.append(cos_sim_iset)
    
# len(cos_sim_sets), len(cos_sim_sets[0]) # n of recording, n of decorrelation pairs (nisi * nori)


# In[863]:


# ## save boot
# # V1_decorr_boot = cos_sim_boots.copy()
# # LM_decorr_boot = cos_sim_boots.copy()
# # LI_decorr_boot = cos_sim_boots.copy()

# # decorr_boot = {'V1_decorr_boot': V1_decorr_boot,
# #                'LM_decorr_boot': LM_decorr_boot,
# #                'LI_decorr_boot': LI_decorr_boot,}

# # import pickle
# # pickle_dir = r'C:\Users\lan\Documents\repos\inter\results\decorrelation vs adp\ref=90 deg'.replace('\\', os.sep)
# # pickle_name = 'decorr_vs_stim2_dist_errorbar_bootstrap_CI_' + '.pickle'
# # with open(os.path.join(pickle_dir, pickle_name), 'wb') as handle:
# #     pickle.dump(decorr_boot, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
# ## load boot
# with open(os.path.join(pickle_dir, pickle_name), 'rb') as handle:
#     decorr_boot = pickle.load(handle)
# decorr_boot.keys()

# V1_decorr_boot = decorr_boot['V1_decorr_boot']
# LM_decorr_boot = decorr_boot['LM_decorr_boot']
# LI_decorr_boot = decorr_boot['LI_decorr_boot']


# In[1054]:


cos_sim_sets = np.array(cos_sim_boots)
cos_sim_sets_avg = np.mean(cos_sim_sets, axis=0) # across boots

perc = 2.5
plt.plot(cos_sim_sets_avg[:nori], color='blue', alpha=0.7)
plt.fill_between(x=np.arange(nori),
                y1=np.percentile(cos_sim_sets, perc, axis=0)[:nori],
                y2=np.percentile(cos_sim_sets, 100-perc, axis=0)[:nori],
                label='isi 250',
                color='blue', alpha=0.2)

plt.plot(cos_sim_sets_avg[nori:], color='orange', alpha=0.7)
plt.fill_between(x=np.arange(nori),
                y1=np.percentile(cos_sim_sets, perc, axis=0)[nori:],
                y2=np.percentile(cos_sim_sets, 100-perc, axis=0)[nori:],
                label='isi 6000',
                color='orange', alpha=0.2);
plt.legend(frameon=False, loc='lower right')

# ## annotate significance
# for iori in np.arange(nori)[:-1]: # exclude 90 deg
#     plt.annotate(sig_star_list[iori], xy=(iori, 1), fontsize=12, ha='center', va='center')

# plt.ylim(0.2, 1.05)
plt.xlabel('stim2 distance from adapter (deg)');
plt.ylabel('cosine similarity')
plt.xticks(range(nori), ['0', '22.5', '45', '67.5', '90'])

fig.tight_layout();
fig_dir = r'C:\Users\lan\Documents\repos\inter\results\decorrelation vs adp\ref=90 deg'.replace('\\', os.sep)
fig_name = f'decorr_fixed_resample_boot_CI{int(100-2*perc)}_median_' + current_area + '.pdf'
plt.savefig(os.path.join(fig_dir, fig_name), bbox_inches='tight')


# ### decorr vs stim2_ori

# In[397]:


df_pop_vec = (df_decorr.groupby(['isi', 'stim2_id', 'cell_id'])
                .dfof.mean().reset_index()
                .pivot_table(index=['isi', 'stim2_id'], columns='cell_id', values='dfof')
                )
df_pop_vec


# In[398]:


nori = df_decorr.stim2_id.nunique()
nisi_now = df_decorr.isi.nunique() # discarded isi 750

cos_sim_list = []
for iisi in range(nisi_now):
    
    if iisi == 0:
        iref = 4 # use 90 deg as ref
    elif iisi == 1:
        iref = nori+4 # -1 (0-based indexing) and +1 (go to next isi) cancel out
    pop_vec_ref = df_pop_vec.iloc[iref, :].values # reference pop vec: 0 deg target resp
    
    for iori in range(nori):
        irow = iisi * nori + iori
        # print(iisi, iori, irow)
        pop_vec_ori = df_pop_vec.iloc[irow, :].values
        cos_sim = dot(pop_vec_ref, pop_vec_ori) / (norm(pop_vec_ref) * norm(pop_vec_ori))
        cos_sim_list.append(cos_sim)


# In[399]:


cos_sim_list = np.array(cos_sim_list)

fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharex=True)
ax[0].plot(cos_sim_list[:nori], color='blue', label='isi 250')
ax[0].plot(cos_sim_list[nori:2*nori], color='orange', label='isi 6000')
ax[0].legend(frameon=False)

ax[1].plot(cos_sim_list[nori:2*nori] - cos_sim_list[:nori], color='red', label='diff')
ax[1].legend(frameon=False)

ax[0].set_xlabel('stim2 distance from adapter (deg)');
ax[0].set_ylabel('cosine similarity')
ax[0].set_xticks(range(nori), ['0', '22.5', '45', '67.5', '90', '112.5', '135', '157.5'], rotation=45)

ax[1].set_xlabel('stim2 distance from adapter (deg)');
ax[1].set_ylabel('cosine similarity')
ax[1].set_xticks(range(nori), ['0', '22.5', '45', '67.5', '90', '112.5', '135', '157.5'], rotation=45)

plt.tight_layout()
fig_dir = r'C:\Users\lan\Documents\repos\inter\results\decorrelation vs adp\ref=90 deg'.replace('\\', os.sep)
fig_name = 'decorr_vs_stim2_ori_' + current_area + '.pdf'
plt.savefig(os.path.join(fig_dir, fig_name), bbox_inches='tight')


# ### decorr vs stim2_distance from adapter

# In[400]:


df_pop_vec = (df_decorr.groupby(['isi', 'stim2_dist', 'cell_id'])
                .dfof.mean().reset_index()
                .pivot_table(index=['isi', 'stim2_dist'], columns='cell_id', values='dfof')
                )
df_pop_vec


# In[401]:


nori = df_decorr.stim2_dist.nunique()
nisi_now = df_decorr.isi.nunique() # discarded isi 750

cos_sim_list = []
for iisi in range(nisi_now):
    
    if iisi == 0:
        iref = 4
    elif iisi == 1:
        iref = nori+4 # -1 (0-based indexing) and +1 (go to next isi) cancel out
    pop_vec_ref = df_pop_vec.iloc[iref, :].values # reference pop vec: 0 deg target resp
    
    for iori in range(nori):
        irow = iisi * nori + iori
        # print(iisi, iori, irow)
        pop_vec_ori = df_pop_vec.iloc[irow, :].values
        cos_sim = dot(pop_vec_ref, pop_vec_ori) / (norm(pop_vec_ref) * norm(pop_vec_ori))
        cos_sim_list.append(cos_sim)


# In[402]:


cos_sim_list = np.array(cos_sim_list)

fig, ax = plt.subplots(1, 2, figsize=(10, 4), sharex=True)
ax[0].plot(cos_sim_list[:nori], color='blue', label='isi 250')
ax[0].plot(cos_sim_list[nori:2*nori], color='orange', label='isi 6000')
ax[0].legend(frameon=False)

ax[1].plot(cos_sim_list[nori:2*nori] - cos_sim_list[:nori], color='red', label='diff')
ax[1].legend(frameon=False)

ax[0].set_xlabel('stim2 distance from adapter (deg)');
ax[0].set_ylabel('cosine similarity')
ax[0].set_xticks(range(nori), ['0', '22.5', '45', '67.5', '90'])

ax[1].set_xlabel('stim2 distance from adapter (deg)');
ax[1].set_ylabel('cosine similarity')
ax[1].set_xticks(range(nori), ['0', '22.5', '45', '67.5', '90'])

plt.tight_layout()
fig_dir = r'C:\Users\lan\Documents\repos\inter\results\decorrelation vs adp\ref=90 deg'.replace('\\', os.sep)
fig_name = 'decorr_vs_stim2_dist_' + current_area + '.pdf'
# plt.savefig(os.path.join(fig_dir, fig_name), bbox_inches='tight')


# # linear SVC decoder
# using function sklearn.svm.LinearSVC:  
# 
# Similar to SVC with parameter kernel=’linear’, but implemented in terms of liblinear rather than libsvm, so it has more flexibility in the choice of penalties and loss functions and should scale better to large numbers of samples.
# 
# This class supports both dense and sparse input and the multiclass support is handled according to a one-vs-the-rest scheme.

# ## example

# In[491]:


from sklearn.svm import LinearSVC # NOTE: current scikit-learn ver: 1.2.2, not newest, beware with doc
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from sklearn.datasets import make_classification

X, y = make_classification(n_features=4, random_state=0)
X.shape, y.shape, y[:5]


# In[492]:


clf = make_pipeline(StandardScaler(),
                    LinearSVC(random_state=0, tol=1e-5))
clf.fit(X, y)

print(clf.named_steps['linearsvc'].coef_)
print(clf.named_steps['linearsvc'].intercept_)

print(clf.predict([[0, 0, 0, 0]]))


# In[490]:


# get accuracy of clf
from sklearn.metrics import accuracy_score

y_pred = clf.predict(X)
acc = accuracy_score(y, y_pred)
print(acc)

# confusion matrix
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay
predictions = clf.predict(X)
cm = confusion_matrix(y, predictions, labels=clf.classes_)
disp = ConfusionMatrixDisplay(confusion_matrix=cm,
                              display_labels=clf.classes_)
disp.plot()

plt.show()


# ## decode 90 vs another ori

# ### filter cell by vis pval

# In[10]:


df_tidy.groupby('cell_id').filter_cell_vis_pval.first().hist(bins=100);
# plt.xlim(0, 0.05);
plt.xlabel('min p-value of cell across stims');
plt.ylabel('cell count');


# In[578]:


pd.options.mode.chained_assignment = None  # default='warn'

df_svc = df_tidy[['dfof', 
                'cell_id', 'resp_id', 'isi', 'stim2_id', 'trial_id',
                'filter_cell_vis', 'filter_cell_stim', 'filter_cell_vis_pval',
                'area', 'mouse', 'date', 'sess'
                ]]

df_svc['stim2_dist'] = df_svc.stim2_id.apply(lambda x: 22.5*(8-x) if x > 4 else 22.5*x) 
                          # stim2 distance from 0 deg, hard coded for 8 ori. TODO: modify if needed
df_svc['date_sess'] = df_svc.date + '_' + df_svc.sess

current_area = 'V1' # TODO: area = subplot
current_datesess = df_svc.date_sess.unique()[0] # TODO: date_sess = errorbar

df_svc = df_svc[(df_svc.isi != 750)
                & (df_svc.filter_cell_vis == True)
                & (df_svc.resp_id == 'R2') # only decode R2, with or without adapter (differentiate by isi)
                & (df_svc.area == current_area)
                & (df_svc.date_sess == current_datesess)
                ]

ncell_keep = 15
vis_pval_thresh = df_svc.groupby('cell_id').first().filter_cell_vis_pval.sort_values()[ncell_keep] # nsample (nrep trial) too small, need to reduce nfeature (ncell)
df_svc = df_svc[df_svc.filter_cell_vis_pval < vis_pval_thresh]

df_svc.cell_id.nunique()


# ### filter cell by SNR

# In[572]:


def calc_SNR(df_tidy):
    ## SNR of R1, aka (R2 without adapter, isi=6k)
    tmp = df_tidy[(df_tidy.resp_id == 'R2') & (df_tidy.isi == 6000)].groupby(['cell_id', 'stim2_id']).dfof
    SNR_R1 = (tmp.mean() / tmp.std()).reset_index()

    ## SNR of R2, aka (R2 with adapter, isi=250)
    tmp = df_tidy[(df_tidy.resp_id == 'R2') & (df_tidy.isi == 250)].groupby(['cell_id', 'stim2_id']).dfof
    SNR_R2 = (tmp.mean() / tmp.std()).reset_index()

    ## merge to df
    df_SNR = SNR_R1.merge(SNR_R2, on=['cell_id', 'stim2_id'], suffixes=('_R1', '_R2'))
    df_SNR.columns = ['cell_id', 'stim2_id', 'SNR_R1', 'SNR_R2'] # rename columns
    df_SNR = df_SNR.sort_values(by='SNR_R1', ascending=False)
    return df_SNR

df_SNR = calc_SNR(df_tidy)


# In[573]:


df_SNR.SNR_R1.hist(bins=100, alpha=0.5);
df_SNR.SNR_R2.hist(bins=100, alpha=0.5);


# In[579]:


df_SNR


# In[574]:


plt.scatter(df_SNR.SNR_R1, df_SNR.SNR_R2, alpha=0.5, s=1);

# correlation
from scipy.stats import pearsonr
pearsonr(df_SNR.SNR_R1, df_SNR.SNR_R2)


# In[575]:


df_SNR.groupby('cell_id').SNR_R2.max().hist(bins=50);
df_SNR.groupby('cell_id').SNR_R2.max().reset_index().sort_values('SNR_R2', ascending=False).head(15).cell_id.values


# In[576]:


pd.options.mode.chained_assignment = None  # default='warn'

df_svc = df_tidy[['dfof', 
                'cell_id', 'resp_id', 'isi', 'stim2_id', 'trial_id',
                'filter_cell_vis', 'filter_cell_stim', 'filter_cell_vis_pval',
                'area', 'mouse', 'date', 'sess'
                ]]

df_svc['stim2_dist'] = df_svc.stim2_id.apply(lambda x: 22.5*(8-x) if x > 4 else 22.5*x) 
                          # stim2 distance from 0 deg, hard coded for 8 ori. TODO: modify if needed
df_svc['date_sess'] = df_svc.date + '_' + df_svc.sess

current_area = 'V1' # TODO: area = subplot
current_datesess = df_svc.date_sess.unique()[0] # TODO: date_sess = errorbar

df_svc = df_svc[(df_svc.isi != 750)
                & (df_svc.filter_cell_vis == True)
                & (df_svc.resp_id == 'R2') # only decode R2, with or without adapter (differentiate by isi)
                & (df_svc.area == current_area)
                & (df_svc.date_sess == current_datesess)
                ]

ncell_keep = 15
df_SNR = calc_SNR(df_svc)
model_cell_id = (df_SNR.groupby('cell_id')
                 .SNR_R2.max().reset_index()
                 .sort_values('SNR_R2', ascending=False)
                 .head(ncell_keep).cell_id.values
                 )

df_svc = df_svc[df_svc.cell_id.isin(model_cell_id)]
assert set(df_svc.cell_id.unique()) == set(sorted(model_cell_id))
df_svc.cell_id.nunique()


# ### validate diversity of tuning curves 
# neuron in model should have diff tuning

# In[577]:


model_cell_resp = df_svc.groupby(['cell_id', 'stim2_id']).dfof.mean().reset_index()
model_cell_std = df_svc.groupby(['cell_id', 'stim2_id']).dfof.std().reset_index()

tuning_pop = []
tuning_pop_std = []
for icell in df_svc.cell_id.unique():
    tuning_cell = model_cell_resp[model_cell_resp.cell_id==icell].dfof.tolist()
    tuning_cell_std = model_cell_std[model_cell_std.cell_id==icell].dfof.tolist()
    
    # tuning_cell = (tuning_cell - np.min(tuning_cell)) / (np.max(tuning_cell) - np.min(tuning_cell)) # min-max normalization
    tuning_pop.append(tuning_cell)
    tuning_pop_std.append(tuning_cell_std)

len(tuning_pop), len(tuning_pop[0]) # ncell_keep, nori
tuning_pop = np.array(tuning_pop)
tuning_pop_std = np.array(tuning_pop_std)
tuning_pop = tuning_pop[np.argsort(np.argmax(tuning_pop, axis=1)), :] # sort cells by tuning argmax
tuning_pop_std = tuning_pop_std[np.argsort(np.argmax(tuning_pop, axis=1)), :]

nsubplot = 5
# fig, ax = plt.subplots(1, nsubplot, figsize=(18, 5), sharex=True, sharey=True)
fig, ax = plt.subplots(1, nsubplot, figsize=(18, 5), sharex=True)
ncell_subplot = ncell_keep // nsubplot

for isubplot in np.arange(nsubplot):
    for icell in np.arange(ncell_subplot*(isubplot-1), ncell_subplot*isubplot):
        ax[isubplot].errorbar(x=np.arange(8) + icell*0.1, 
                              y=tuning_pop[icell], 
                              yerr=tuning_pop_std[icell], 
                              alpha=0.5)
        
fig.text(0.5, -0.05, 'stim2 orientation', ha='center', fontsize=18)
fig.text(-0.02, 0.5, 'normed tuning', va='center', rotation='vertical', fontsize=18)
# fig.text(-0.02, 0.5, 'tuning', va='center', rotation='vertical', fontsize=18)

fig.tight_layout()
# fig.savefig(dir_fig + 'tuning_of_model_cells_no_norm.pdf', bbox_inches='tight')


# In[402]:


def df_to_train_test(df):
    # input: filtered df_svc
    # output: X_train, X_test, y_train, y_test (of the filter condition)
    
    label_arr = df.groupby('trial_id').stim2_dist.first().values
    label_arr = (label_arr == 90) # 1 = 90 deg, 0 = other ori
    
    feature_mat = (df
                   .pivot_table(index=['trial_id'], columns='cell_id', values='dfof')
                   .fillna(0).to_numpy())
    
    assert feature_mat.shape[0] == label_arr.shape[0]
    # print(label_arr.shape) # ntrial
    # print(feature_mat.shape) # ntrial x ncell
    
    X_train, X_test, y_train, y_test = train_test_split(
        feature_mat, label_arr, test_size=int(1), # 1 trial for test. do not differentiate val/test now
        random_state=0, shuffle=True) 
    # print(X_train.shape, y_train.shape, X_test.shape, y_test.shape) # ntrial x ncell
    
    return X_train, X_test, y_train, y_test


# ### try

# In[88]:


acc_train_250_arr = []
acc_train_6000_arr = []
acc_test_250_arr = []
acc_test_6000_arr = []

# other_dist = sorted(df_svc.stim2_dist.unique())[:-1][::-1] # [0, 22, 45, 67] -> [67.5, 45.0, 22.5, 0.0]
tmp = df_svc.stim2_id.unique()
other_ori = sorted(tmp[tmp != 4])

# for idist in other_dist: # exclude 90 deg
for iori in other_ori: # exclude 90 deg
    df_pair = df_svc[((df_svc.stim2_id == 4) | (df_svc.stim2_id == iori))] # 90 deg vs another ori
    
    df_pair_250 = df_pair[df_pair.isi == 250]
    df_pair_6000 = df_pair[df_pair.isi == 6000]
    
    ## split train vs test for each isi
    X_train_250, X_test_250, y_train_250, y_test_250 = df_to_train_test(df_pair_250)
    X_train_6000, X_test_6000, y_train_6000, y_test_6000 = df_to_train_test(df_pair_6000)
    
    ## stack training data across isi - ntrial per isi mostly balanced (30% vs 35%)
    X_train = np.vstack((X_train_250, X_train_6000))
    y_train = np.hstack((y_train_250, y_train_6000))
    
    ## shuffle trials in the same way
    np.random.seed(42)
    np.random.shuffle(X_train) # shuffled along the first axis, aka trials
    np.random.shuffle(y_train)
    
    ## fit train (merge isi)
    clf = make_pipeline(StandardScaler(),
                        LinearSVC(random_state=42, tol=1e-5, max_iter=100000, 
                                  penalty='l1', dual=False, C=0.1))
    clf.fit(X_train, y_train)
    
    ## predict test (separate isi)
    # y_pred_train_250 = clf.predict(X_train_250)
    # y_pred_train_6000 = clf.predict(X_train_6000)
    # y_pred_test_250 = clf.predict(X_test_250)
    # y_pred_test_6000 = clf.predict(X_test_6000)
    
    acc_train_250 = clf.score(X_train_250, y_train_250)
    acc_train_6000 = clf.score(X_train_6000, y_train_6000)
    acc_test_250 = clf.score(X_test_250, y_test_250)
    acc_test_6000 = clf.score(X_test_6000, y_test_6000)
    
    acc_train_250_arr.append(acc_train_250)
    acc_train_6000_arr.append(acc_train_6000)
    acc_test_250_arr.append(acc_test_250)
    acc_test_6000_arr.append(acc_test_6000)
    
    # ## confusion matrix
    # cm_train = confusion_matrix(y_train, y_pred_train, labels=clf.classes_)
    # cm_test_250 = confusion_matrix(y_test_250, y_pred_250, labels=clf.classes_)
    # cm_test_6000 = confusion_matrix(y_test_6000, y_pred_6000, labels=clf.classes_)
    
    # break
    
X_train.shape, y_train.shape, #X_test_250.shape, y_test_250.shape, X_test_6000.shape, y_test_6000.shape
# disp = ConfusionMatrixDisplay(confusion_matrix=cm_test_6000,
#                               display_labels=clf.classes_)
# disp.plot()
# plt.show()


# In[89]:


x = other_ori
plt.plot(x, acc_train_250_arr, label='train 250')
plt.plot(x, acc_train_6000_arr, label='train 6000')
plt.plot(x, acc_test_250_arr, label='test 250')
plt.plot(x, acc_test_6000_arr, label='test 6000')

plt.legend(frameon=False);
plt.xticks(x);
plt.xlabel('stim2 id');
plt.ylabel('accuracy');


# In[36]:


from sklearn import datasets
from sklearn.model_selection import cross_val_score
from sklearn import svm

# X, y = datasets.load_iris(return_X_y=True)
# clf = svm.SVC(kernel='linear', C=1, random_state=42)
scores = cross_val_score(clf, X_train, y_train, cv=10)
print(scores)
clf.set_params(linearsvc__C=0.001)
clf.get_params()['linearsvc__C']


# In[38]:


from sklearn.model_selection import cross_validate
# from sklearn.metrics import recall_score

# scoring = ['precision_macro', 'recall_macro']
clf = svm.SVC(kernel='linear', C=1, random_state=0)
scores = cross_validate(clf, X_train, y_train, cv=5,
                        # scoring=scoring, 
                        return_train_score=True)
scores['test_score'], scores['train_score']


# In[239]:


from sklearn.model_selection import LeaveOneOut

X = [1, 2, 3, 4]
loo = LeaveOneOut()
for train, test in loo.split(X):
    print("%s %s" % (train, test))


# ## cross val to search reg param

# ### aware decoder
# 
# two separate models for two isi.  
# one model: train and test with isi=6k  
# another model: train and test with isi=250  

# In[560]:


## find optimal reg param C. 

df_svc = get_df_svc(df_tidy, iarea='V1', idatesess=4, ncell_keep=999)

# use_data = 'isi_both'
# use_data = 'isi=250' # adapted
use_data = 'isi=6k' # unadapted

nrow = 4
ncol = 2
cross_val_method = LeaveOneOut()
fig, ax = plt.subplots(4, 2, figsize=(12, 20), sharey=True)

acc_val_mean_iori = []
acc_val_ci_iori = [] # confidence interval based on binomial distribution
for iori in tqdm(other_ori):
    df_pair = df_svc[((df_svc.stim2_id == 4) | (df_svc.stim2_id == iori))]
    df_pair_250 = df_pair[df_pair.isi == 250]
    df_pair_6000 = df_pair[df_pair.isi == 6000]

    ## split train vs test for each isi
    X_train_250, X_test_250, y_train_250, y_test_250 = df_to_train_test(df_pair_250)
    X_train_6000, X_test_6000, y_train_6000, y_test_6000 = df_to_train_test(df_pair_6000)

    ## stack training data across isi - ntrial per isi mostly balanced (30% vs 35%)
    X_train = np.vstack((X_train_250, X_train_6000))
    y_train = np.hstack((y_train_250, y_train_6000))
    X_test = np.vstack((X_test_250, X_test_6000))
    y_test = np.hstack((y_test_250, y_test_6000))

    ## shuffle trials in the same way
    seed_val = 0
    np.random.seed(seed_val)
    np.random.shuffle(X_train) # shuffled along the first axis, aka trials
    np.random.shuffle(y_train)
    np.random.shuffle(X_test)
    np.random.shuffle(y_test)

    ## cross validation to find optimal C
    acc_train_mean = []
    acc_train_std = []
    acc_val_mean = []
    acc_val_std = []
    acc_val_ci = [] # confidence interval based on binomial distribution
    acc_test = []
    # C_list = [1e-4 * 10**i for i in range(8)] # 1e-4, 1e-3, ..., 1e3, 1e4
    C_list = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1] # zoom in around 0.1
    nfold = 5

    for C_val in C_list:
        clf = make_pipeline(StandardScaler(),
                            LinearSVC(tol=1e-5, max_iter=int(1e6), 
                                    penalty='l1', dual=False, C=C_val))
        
        if use_data == 'isi_both':
            score_val = cross_validate(clf, X_train, y_train, cv=cross_val_method, return_train_score=True) # instead of cv=nfold, we use leave-one-out
            clf.fit(X_train, y_train)
            score_test = clf.score(X_test, y_test)
        elif use_data == 'isi=250':
            score_val = cross_validate(clf, X_train_250, y_train_250, cv=cross_val_method, return_train_score=True)
            clf.fit(X_train_250, y_train_250)
            score_test = clf.score(X_test_250, y_test_250)
        elif use_data == 'isi=6k':
            score_val = cross_validate(clf, X_train_6000, y_train_6000, cv=cross_val_method, return_train_score=True)
            clf.fit(X_train_6000, y_train_6000)
            score_test = clf.score(X_test_6000, y_test_6000)
        
        acc_train_mean.append(np.mean(score_val['train_score']))
        acc_train_std.append(np.std(score_val['train_score']))
        acc_val_mean.append(np.mean(score_val['test_score']))
        acc_val_std.append(np.std(score_val['test_score']))
        ci_low, ci_high = proportion_confint(count=score_val['test_score'].sum(), nobs=len(score_val['test_score']))
        ci_err = (ci_high - ci_low)/2
        acc_val_ci.append(ci_err)
        acc_test.append(score_test)
    acc_val_mean_iori.append(acc_val_mean)
    acc_val_ci_iori.append(acc_val_ci)
        
        
    ax[int(iori/ncol), iori%ncol].errorbar(x=C_list, y=acc_train_mean, 
                yerr=acc_train_std, fmt='o', capsize=5, capthick=2, label='train')
    ax[int(iori/ncol), iori%ncol].errorbar(x=C_list, y=acc_val_mean,
                yerr=acc_val_ci, fmt='o', capsize=5, capthick=2, label='val')
    # ax[int(iori/ncol), iori%ncol].plot(C_list, acc_test, color='r', marker=".", alpha=0.5, label='test')
    ax[int(iori/ncol), iori%ncol].axhline(y=0.5, color='k',alpha=0.2, label='chance')
    ax[int(iori/ncol), iori%ncol].axhline(y=1, color='k',alpha=0.2, label='1')

    # set x axis to log scale
    ax[int(iori/ncol), iori%ncol].set_xscale('log')

    ax[int(iori/ncol), iori%ncol].legend(frameon=False);
    ax[int(iori/ncol), iori%ncol].set_title(f'decoder: {int(22.5*iori)} vs 90 deg');
    ax[int(iori/ncol), iori%ncol].set_xlabel('more reg ... <-- C --> ... less reg');
    ax[int(iori/ncol), iori%ncol].set_ylabel('accuracy');
    
fig.tight_layout()
print(use_data)

dir_fig = r'C:\Users\lan\Documents\repos\inter\results\decoder_grat8/'.replace('\\', '/')
# fig.savefig(dir_fig + f'aware_decoder_L1_{use_data}_all_cell.pdf', bbox_inches='tight')


# In[564]:


fix_C = 0.1
# if use_data == 'isi=6k':
    # fix_C = 0.1
# if use_data == 'isi=250':
    # fix_C = 0.2
    
acc_val_C = np.array(acc_val_mean_iori)[:, np.where(np.array(C_list)==fix_C)[0][0]]
acc_val_ci_C = np.array(acc_val_ci_iori)[:, np.where(np.array(C_list)==fix_C)[0][0]]

plt.errorbar(x=np.arange(len(stim2_ori_list)), 
             y=acc_val_C, 
             yerr=acc_val_ci_C, 
             fmt='o', color='b', alpha=0.5, ecolor='lightgray', elinewidth=3, capsize=0);

stim2_ori_list = ['0', '22.5', '45', '67.5', '112.5', '135', '157.5']
plt.xticks(np.arange(len(stim2_ori_list)), stim2_ori_list);
plt.xlabel('stim2_ori')
plt.ylabel(f'acc_val');

title_str = f'aware decoder, leave one out, binofit, {use_data}, C={fix_C}'
plt.title(title_str);
# plt.savefig(dir_fig + f'{title_str}.pdf', bbox_inches='tight')


# ### naive unaware decoder
# train with 6k, test with either - more directly related to decorr, give lower bound of decoding acc

# In[289]:


## find optimal reg param C. 

nrow = 4
ncol = 2
fig, ax = plt.subplots(4, 2, figsize=(12, 20), sharey=True)

acc_val_mean_iori = []
acc_val_ci_iori = [] # confidence interval based on binomial distribution
acc_test_iori = []

for iori in tqdm(other_ori):
    
    ## get data for each pair of ori (another vs 90 deg)
    df_pair = df_svc[((df_svc.stim2_id == 4) | (df_svc.stim2_id == iori))]
    df_pair_250 = df_pair[df_pair.isi == 250]
    df_pair_6000 = df_pair[df_pair.isi == 6000]

    ## split train vs test for each isi
    X_train_250, X_test_250, y_train_250, y_test_250 = df_to_train_test(df_pair_250)
    X_train_6000, X_test_6000, y_train_6000, y_test_6000 = df_to_train_test(df_pair_6000)

    # ## stack training data across isi - ntrial per isi mostly balanced (30% vs 35%)
    # X_train = np.vstack((X_train_250, X_train_6000))
    # y_train = np.hstack((y_train_250, y_train_6000))
    # X_test = np.vstack((X_test_250, X_test_6000))
    # y_test = np.hstack((y_test_250, y_test_6000))

    # ## shuffle trials in the same way
    # seed_val = 0
    # np.random.seed(seed_val)
    # np.random.shuffle(X_train) # shuffled along the first axis, aka trials
    # np.random.shuffle(y_train)
    # np.random.shuffle(X_test)
    # np.random.shuffle(y_test)

    ## cross validation to find optimal C
    acc_train_mean = []
    acc_train_std = []
    acc_val_mean = []
    acc_val_std = []
    acc_val_ci = [] # confidence interval based on binomial distribution
    acc_test = []
    # C_list = [1e-4 * 10**i for i in range(8)] # 1e-4, 1e-3, ..., 1e3, 1e4
    C_list = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1] # zoom in around 0.1
    nfold = 5

    # use_data = 'isi=250' # train with unadapted, test with adapted
    # use_data = 'isi=6k' # train with unadapted, test with unadapted


    for C_val in C_list:
        clf = make_pipeline(StandardScaler(),
                            LinearSVC(tol=1e-5, max_iter=int(1e6), 
                                    penalty='l1', dual=False, C=C_val))

        ## train and val with isi 6k. val is considered as test
        score_val = cross_validate(clf, X_train_6000, y_train_6000, cv=LeaveOneOut(), return_train_score=True)        
        acc_train_mean.append(np.mean(score_val['train_score']))
        acc_train_std.append(np.std(score_val['train_score']))
        
        acc_val_mean.append(np.mean(score_val['test_score']))
        acc_val_std.append(np.std(score_val['test_score']))
        ci_low, ci_high = proportion_confint(count=score_val['test_score'].sum(), nobs=len(score_val['test_score']))
        ci_err = (ci_high - ci_low)/2
        acc_val_ci.append(ci_err)

        ## test with isi 250. train_250 is considered as test (bc model is trained with isi 6k)
        clf.fit(X_train_6000, y_train_6000)
        score_test = clf.score(X_train_250, y_train_250) # test data has only 1 trial. use train as test for isi 250
        acc_test.append(score_test)
        
    acc_val_mean_iori.append(acc_val_mean)
    acc_val_ci_iori.append(acc_val_ci)
    acc_test_iori.append(acc_test)
        
        
    ax[int(iori/ncol), iori%ncol].errorbar(x=C_list, y=acc_train_mean, 
                yerr=acc_train_std, fmt='o', capsize=5, capthick=2, label='train 6k')
    ax[int(iori/ncol), iori%ncol].errorbar(x=C_list, y=acc_val_mean,
                yerr=acc_val_ci, fmt='o', capsize=5, capthick=2, label='val 6k')
    ax[int(iori/ncol), iori%ncol].plot(C_list, acc_test, color='r', marker=".", alpha=0.5, label='test 250')
    ax[int(iori/ncol), iori%ncol].axhline(y=0.5, color='k',alpha=0.2, label='chance')

    # set x axis to log scale
    ax[int(iori/ncol), iori%ncol].set_xscale('log')

    ax[int(iori/ncol), iori%ncol].legend(frameon=False);
    ax[int(iori/ncol), iori%ncol].set_title(f'decoder: {int(22.5*iori)} vs 90 deg');
    ax[int(iori/ncol), iori%ncol].set_xlabel('more reg ... <-- C --> ... less reg');
    ax[int(iori/ncol), iori%ncol].set_ylabel('accuracy');
    
fig.tight_layout()

dir_fig = r'C:\Users\lan\Documents\repos\inter\results\decoder_grat8/'.replace('\\', '/')
fig.savefig(dir_fig + f'decoder_L1_train_val_6k_test_250.pdf', bbox_inches='tight')


# In[294]:


fix_C = 0.1
# if use_data == 'isi=6k':
    # fix_C = 0.1
# if use_data == 'isi=250':
    # fix_C = 0.2
    
acc_val_C = np.array(acc_val_mean_iori)[:, np.where(np.array(C_list)==fix_C)[0][0]]
acc_val_ci_C = np.array(acc_val_ci_iori)[:, np.where(np.array(C_list)==fix_C)[0][0]]
acc_test_C = np.array(acc_test_iori)[:, np.where(np.array(C_list)==fix_C)[0][0]]

plt.errorbar(x=np.arange(len(stim2_ori_list)), 
             y=acc_val_C, 
             yerr=acc_val_ci_C, 
             label='acc val isi=6k',
             fmt='o', color='b', alpha=0.5, ecolor='lightgray', elinewidth=3, capsize=0);
plt.plot(acc_test_C, 'o', color='r', alpha=0.5, label='acc test isi=250')

stim2_ori_list = ['0', '22.5', '45', '67.5', '112.5', '135', '157.5']
plt.xticks(np.arange(len(stim2_ori_list)), stim2_ori_list);
plt.xlabel('stim2_ori')
plt.ylabel(f'acc');
plt.legend(frameon=False, loc='lower left')

title_str = f'naive unaware decoder, binofit, C={fix_C}'
plt.title(title_str);
plt.savefig(dir_fig + f'{title_str}.pdf', bbox_inches='tight')


# ### experienced unaware decoder
# train with both, test with either - assume decoder learned from visual xp

# In[297]:


## find optimal reg param C. 

nrow = 4
ncol = 2
cross_val_method = LeaveOneOut()
fig, ax = plt.subplots(4, 2, figsize=(12, 20), sharey=True)

acc_val_mean_iori = []
acc_val_ci_iori = [] # confidence interval based on binomial distribution
for iori in tqdm(other_ori):
    df_pair = df_svc[((df_svc.stim2_id == 4) | (df_svc.stim2_id == iori))]
    df_pair_250 = df_pair[df_pair.isi == 250]
    df_pair_6000 = df_pair[df_pair.isi == 6000]

    ## split train vs test for each isi
    X_train_250, X_test_250, y_train_250, y_test_250 = df_to_train_test(df_pair_250)
    X_train_6000, X_test_6000, y_train_6000, y_test_6000 = df_to_train_test(df_pair_6000)

    ## stack training data across isi - ntrial per isi mostly balanced (30% vs 35%)
    X_train = np.vstack((X_train_250, X_train_6000))
    y_train = np.hstack((y_train_250, y_train_6000))
    X_test = np.vstack((X_test_250, X_test_6000))
    y_test = np.hstack((y_test_250, y_test_6000))

    # ## shuffle trials in the same way
    # seed_val = 0
    # np.random.seed(seed_val)
    # np.random.shuffle(X_train) # shuffled along the first axis, aka trials
    # np.random.shuffle(y_train)
    # np.random.shuffle(X_test)
    # np.random.shuffle(y_test)

    ## cross validation to find optimal C
    acc_train_mean = []
    acc_train_std = []
    acc_val_mean = []
    acc_val_std = []
    acc_val_ci = [] # confidence interval based on binomial distribution
    acc_test = []
    # C_list = [1e-4 * 10**i for i in range(8)] # 1e-4, 1e-3, ..., 1e3, 1e4
    C_list = [0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1] # zoom in around 0.1
    nfold = 5

    use_data = 'isi_both'
    # use_data = 'isi=250' # adapted
    # use_data = 'isi=6k' # unadapted


    for C_val in C_list:
        clf = make_pipeline(StandardScaler(),
                            LinearSVC(tol=1e-5, max_iter=int(1e6), 
                                    penalty='l1', dual=False, C=C_val))
        
        if use_data == 'isi_both':
            score_val = cross_validate(clf, X_train, y_train, cv=cross_val_method, return_train_score=True) # instead of cv=nfold, use leave-one-out
        #     clf.fit(X_train, y_train)
        #     score_test = clf.score(X_test, y_test)
        # elif use_data == 'isi=250':
        #     score_val = cross_validate(clf, X_train_250, y_train_250, cv=cross_val_method, return_train_score=True)
        #     clf.fit(X_train_250, y_train_250)
        #     score_test = clf.score(X_test_250, y_test_250)
        # elif use_data == 'isi=6k':
        #     score_val = cross_validate(clf, X_train_6000, y_train_6000, cv=cross_val_method, return_train_score=True)
        #     clf.fit(X_train_6000, y_train_6000)
        #     score_test = clf.score(X_test_6000, y_test_6000)
        
        acc_train_mean.append(np.mean(score_val['train_score']))
        acc_train_std.append(np.std(score_val['train_score']))
        acc_val_mean.append(np.mean(score_val['test_score']))
        acc_val_std.append(np.std(score_val['test_score']))
        ci_low, ci_high = proportion_confint(count=score_val['test_score'].sum(), nobs=len(score_val['test_score']))
        ci_err = (ci_high - ci_low)/2
        acc_val_ci.append(ci_err)
        # acc_test.append(score_test)
    acc_val_mean_iori.append(acc_val_mean)
    acc_val_ci_iori.append(acc_val_ci)
        
        
    ax[int(iori/ncol), iori%ncol].errorbar(x=C_list, y=acc_train_mean, 
                yerr=acc_train_std, fmt='o', capsize=5, capthick=2, label='train')
    ax[int(iori/ncol), iori%ncol].errorbar(x=C_list, y=acc_val_mean,
                yerr=acc_val_ci, fmt='o', capsize=5, capthick=2, label='val')
    # ax[int(iori/ncol), iori%ncol].plot(C_list, acc_test, color='r', marker=".", alpha=0.5, label='test')
    ax[int(iori/ncol), iori%ncol].axhline(y=0.5, color='k',alpha=0.2, label='chance')

    # set x axis to log scale
    ax[int(iori/ncol), iori%ncol].set_xscale('log')

    ax[int(iori/ncol), iori%ncol].legend(frameon=False);
    ax[int(iori/ncol), iori%ncol].set_title(f'decoder: {int(22.5*iori)} vs 90 deg');
    ax[int(iori/ncol), iori%ncol].set_xlabel('more reg ... <-- C --> ... less reg');
    ax[int(iori/ncol), iori%ncol].set_ylabel('accuracy');
    
fig.tight_layout()

dir_fig = r'C:\Users\lan\Documents\repos\inter\results\decoder_grat8/'.replace('\\', '/')
fig.savefig(dir_fig + f'decoder_L1_leave1out.pdf', bbox_inches='tight')


# In[298]:


fix_C = 0.1
    
acc_val_C = np.array(acc_val_mean_iori)[:, np.where(np.array(C_list)==fix_C)[0][0]]
acc_val_ci_C = np.array(acc_val_ci_iori)[:, np.where(np.array(C_list)==fix_C)[0][0]]

plt.errorbar(x=np.arange(len(stim2_ori_list)), 
             y=acc_val_C, 
             yerr=acc_val_ci_C, 
             fmt='o', color='b', alpha=0.5, ecolor='lightgray', elinewidth=3, capsize=0);

stim2_ori_list = ['0', '22.5', '45', '67.5', '112.5', '135', '157.5']
plt.xticks(np.arange(len(stim2_ori_list)), stim2_ori_list);
plt.xlabel('stim2_ori')
plt.ylabel(f'acc_val');

title_str = f'experienced unaware decoder, leave one out, binofit, C={fix_C}'
plt.title(title_str);
plt.savefig(dir_fig + f'{title_str}.pdf', bbox_inches='tight')


# # expand
# area subplot, isi color, date_sess errorbar  
# fix C=0.1

# In[504]:


# pd.options.mode.chained_assignment = None  # default='warn'

def get_df_svc(df_tidy, iarea, idatesess, ncell_keep=15):
    df_svc = df_tidy[['dfof', 
                    'cell_id', 'resp_id', 'isi', 'stim2_id', 'trial_id',
                    'filter_cell_vis', 'filter_cell_stim', 'filter_cell_vis_pval',
                    'area', 'mouse', 'date', 'sess'
                    ]]

    df_svc['stim2_dist'] = df_svc.stim2_id.apply(lambda x: 22.5*(8-x) if x > 4 else 22.5*x) 
                            # stim2 distance from 0 deg, hard coded for 8 ori. TODO: modify if needed
    df_svc['date_sess'] = df_svc.date + '_' + df_svc.sess
    
    df_svc = df_svc[(df_svc.isi != 750)
                    & (df_svc.filter_cell_vis == True)
                    & (df_svc.resp_id == 'R2') # only decode R2, with or without adapter (differentiate by isi)
                    & (df_svc.area == iarea)
                    ]
    date_sess_now = df_svc.date_sess.unique()[idatesess] # sess id in a specific area
    df_svc = df_svc[df_svc.date_sess == date_sess_now]

    if df_svc.cell_id.nunique() > ncell_keep:
        vis_pval_thresh = (df_svc.groupby('cell_id').first()
                        .filter_cell_vis_pval.sort_values() # ascending, smaller pval better
                        [ncell_keep]) # nsample (nrep trial) too small, need to reduce nfeature (ncell)
        df_svc = df_svc[df_svc.filter_cell_vis_pval < vis_pval_thresh]
        
    return df_svc


# In[556]:


iarea = 'LI'

true_pos_isi = []
n_sample_isi = []
for use_data in ['isi=250', 'isi=6k']:

    true_pos_ori = np.zeros(len(other_ori) + 1) # true positive across ori pairs. leave a blank for stim2_id=4 (90 deg)
    n_sample_ori = np.zeros(len(other_ori) + 1)

    nsess_iarea = (df_tidy[(df_tidy.area == iarea) 
                           & (df_tidy.filter_cell_vis == True)] # some LI sessions have no vis cell
                .groupby(['date', 'sess'])
                .first().__len__())

    for isess in tqdm(range(nsess_iarea)):
        del df_svc
        df_svc = get_df_svc(df_tidy, iarea=iarea, idatesess=isess, ncell_keep=999)

        for iori in (other_ori):
            
            df_pair = df_svc[((df_svc.stim2_id == 4) | (df_svc.stim2_id == iori))]
            df_pair_250 = df_pair[df_pair.isi == 250]
            df_pair_6000 = df_pair[df_pair.isi == 6000]

            X_train_250, X_test_250, y_train_250, y_test_250 = df_to_train_test(df_pair_250)
            X_train_6000, X_test_6000, y_train_6000, y_test_6000 = df_to_train_test(df_pair_6000)


            clf = make_pipeline(StandardScaler(),
                                LinearSVC(tol=1e-5, max_iter=int(1e6), 
                                        penalty='l1', dual=False, C=0.1))

            if use_data == 'isi=250':
                score_val = cross_validate(clf, X_train_250, y_train_250, cv=LeaveOneOut(), return_train_score=True)
            elif use_data == 'isi=6k':
                score_val = cross_validate(clf, X_train_6000, y_train_6000, cv=LeaveOneOut(), return_train_score=True)

            true_pos = score_val['test_score'].sum()
            n_sample = len(score_val['test_score'])
            true_pos_ori[iori] += true_pos
            n_sample_ori[iori] += n_sample
    
    true_pos_isi.append(true_pos_ori)
    n_sample_isi.append(n_sample_ori)


# In[557]:


# decode_res_V1 = [true_pos_isi, n_sample_isi]
# decode_res_LM = [true_pos_isi, n_sample_isi]
# decode_res_LI = [true_pos_isi, n_sample_isi]

# grat8_decode_all_cells = {'decode_res_V1': decode_res_V1,
#                     'decode_res_LM': decode_res_LM,
#                     'decode_res_LI': decode_res_LI}
# import pickle
# with open(dir_fig + 'grat8_decode_all_cells.pickle', 'wb') as f:
#     pickle.dump(grat8_decode_all_cells, f)


# iarea = 'LM'
# true_pos_isi, n_sample_isi = decode_res_LM


# ## plot

# In[558]:


# suppress RuntimeWarning: invalid value encountered in true_divide
import warnings
warnings.filterwarnings('ignore')

isi_list = ['isi=250', 'isi=6k']
nisi = len(isi_list)
for iisi in range(nisi):

    ci_low, ci_high = proportion_confint(count=true_pos_isi[iisi], 
                                         nobs=n_sample_isi[iisi])
    ci_err = (ci_high - ci_low)/2
    
    plt.errorbar(x=np.arange(8) + iisi*0.1, 
                y=true_pos_isi[iisi] / n_sample_isi[iisi], 
                yerr=ci_err, 
                label=isi_list[iisi], fmt='o', alpha=0.5);

plt.axvline(x=4, color='gray', linestyle='-');
plt.ylim(0, 1);
plt.legend(frameon=False, loc='lower left');
plt.title(iarea);

