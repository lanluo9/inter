#!/usr/bin/env python
# coding: utf-8

# In[4]:


import numpy as np
import pandas as pd
import scipy.io
from scipy import stats

# import seaborn as sns
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

from tqdm import tqdm
import os
from scipy.io import savemat
import pickle
import time
from IPython.display import clear_output

import logging
from datetime import datetime


# # cell filter
# visually responsive & image driven 

# In[5]:


dir_inter = r'Z:\All_Staff\home\lan\Data\2P_images\mat_inter/'.replace('\\', '/')
dir_file = dir_inter + 'adp_dataset_master.xlsx'
data_info = pd.read_excel(dir_file)
data_info.tail()


# In[6]:


meta = data_info[(data_info.paradigm == 'grating') # grating_8ori_multisess
                 & (data_info.gcamp == '6s') # avoid mixing in gcamp8f
                 & ((data_info.cellpose_seg == True) | (data_info.manual_seg == True)) # ensure segmentation
                 ]

meta = meta.reset_index(drop=True)
nset = meta.shape[0]
print(meta.area.value_counts())
# meta


# In[7]:


session_mode = 'single session' # then irun must be appended to dir_sub !
# session_mode = 'multi session' # no need to append irun to dir_sub

## set up logging
plt.set_loglevel(level='warning') # turn off matplotlib debug logging
logging.basicConfig(filename=r'C:\Users\ll357\Documents\inter\data\vis_driven.log'.replace('\\', '/'), level=logging.DEBUG)
logging.info('\n\n\n\n\n') # log several linebreaks to separate different runs
logging.info(str(datetime.now()))

for i in np.arange(nset): # TODO: this overwrites vis-driven.pickle multiple times for multisess data, fix later
    # clear_output(wait=True)
    logging.info(f'set {i+1} of {nset}')
    mouse = meta.iloc[i].mouse.astype(str)
    date = meta.iloc[i].date.astype(str)
    area = meta.iloc[i].area
    irun = meta.iloc[i].num # already a string due to concat lindsey's grating 8ori data
    logging.info(f'{mouse} {date} {irun} {area}')

    ## check if resp_base_trialwise exists
    ## load data
    logging.info('waiting for mat')
    print('waiting for mat')
    print('------------------')

    while True:
        try:
            dir_sub = area + '_i' + mouse + '_' + date
            if session_mode == 'single session':
                dir_sub += '_' + irun # else: multi sess, dont append irun
            print(dir_sub)
            try: # manual seg data has no suffix '_cellpose'
                dfof_trialwise = scipy.io.loadmat(os.path.join(dir_inter, dir_sub, 'resp_base_trialwise' + '.mat'))
            except: # cellpose seg data has suffix '_cellpose'
                dir_sub += '_cellpose'
                dfof_trialwise = scipy.io.loadmat(os.path.join(dir_inter, dir_sub, 'resp_base_trialwise' + '.mat'))
            logging.info('mat loaded')
            break
        except:
            time.sleep(60)
            clear_output(wait=True)
            continue

    if len(dfof_trialwise['dfof_ad_trial'].shape) == 2: # if contains only one ISI, likely 250
        isi_mode = 'only one isi'
        print(f'isi_mode is {isi_mode} (only one isi, grat_SF or bunny)')
        print(dfof_trialwise['dfof_ad_trial'].shape, dfof_trialwise['dfof_ad_trial'][0,0].shape)
        dfof_ad_trial = dfof_trialwise['dfof_ad_trial'] # do not subtract baseline here: stim resp will be compared to base resp
        dfof_base_trial = dfof_trialwise['dfof_base_trial']

    if len(dfof_trialwise['dfof_ad_trial'].shape) == 3: # if contains multiple ISI, likely inf-750-250
        isi_mode = 'grat_8ori_3isi'
        print(f'isi_mode is {isi_mode} (grating with 8 orientations and 3 ISI)')
        print(dfof_trialwise['dfof_ad_trial'].shape, dfof_trialwise['dfof_ad_trial'][0,0,0].shape)
        dfof_tg_trial = dfof_trialwise['dfof_tg_trial'][:,:,0] # only use the first ISI, aka no adapter trials
        dfof_base2_trial = dfof_trialwise['dfof_base2_trial'][:,:,0]
        print('to find visually driven cells and image driven cells, we should take only no-adapter trials and use target response. \n 8 possible target orientations will cover all needed cells')

    ## stats
    if isi_mode == 'only one isi':
        print(isi_mode)

        ncell = dfof_ad_trial.shape[0]
        nstim = dfof_ad_trial.shape[1]

        ## according to Ohki 2020
        # p_anova = np.ones((ncell, 1)) * np.nan # anova -> visually driven cells,  
        # p_kruskal = np.ones((ncell, 1)) * np.nan # t test -> certain image responsive cells,
        p_ttest = np.ones((ncell, nstim)) * np.nan
        evoked = np.ones((ncell, nstim)) * np.nan # amp thresh -> lower false positive rate

        for icell in tqdm(np.arange(ncell)):
            base_cell_anova = np.concatenate([dfof_base_trial[icell, stim] for stim in range(nstim)])
            stim_cell_anova = [] 
            for istim in np.arange(nstim):
                stim_cell_anova.append(np.array(dfof_ad_trial[icell, istim]).flatten())

                base_cell = dfof_base_trial[icell, istim]
                stim_cell = dfof_ad_trial[icell, istim]
                _, p_ttest_i = stats.ttest_ind(base_cell.flatten(), stim_cell.flatten(), equal_var=False, alternative='less')
                p_ttest[icell, istim] = p_ttest_i

                evoked_cell = (stim_cell - base_cell) / (base_cell + 1e-7)
                evoked_i = np.mean(evoked_cell, axis=0) # trial averaged evoked resp
                evoked[icell, istim] = evoked_i
            
            # _, p_anova_cell = stats.f_oneway(np.array(base_cell_anova).flatten(), *stim_cell_anova)
            # p_anova[icell] = p_anova_cell
            # _, p_kruskal_cell = stats.kruskal(np.array(base_cell_anova).flatten(), *stim_cell_anova)
            # p_kruskal[icell] = p_kruskal_cell

    elif isi_mode == 'grat_8ori_3isi':
        print(isi_mode)

        ncell = dfof_tg_trial.shape[0]
        nstim = dfof_tg_trial.shape[1]

        # according to Ohki 2020
        # p_anova = np.ones((ncell, 1)) * np.nan # anova -> visually driven cells,  
        # p_kruskal = np.ones((ncell, 1)) * np.nan # t test -> certain image responsive cells,
        p_ttest = np.ones((ncell, nstim)) * np.nan
        evoked = np.ones((ncell, nstim)) * np.nan # amp thresh -> lower false positive rate

        for icell in tqdm(np.arange(ncell)):
            base_cell_anova = np.concatenate([dfof_base2_trial[icell, stim] for stim in range(nstim)])
            stim_cell_anova = [] 
            for istim in np.arange(nstim):
                stim_cell_anova.append(np.array(dfof_tg_trial[icell, istim]).flatten())

                base_cell = dfof_base2_trial[icell, istim]
                stim_cell = dfof_tg_trial[icell, istim]
                _, p_ttest_i = stats.ttest_ind(base_cell.flatten(), stim_cell.flatten(), equal_var=False, alternative='less')
                p_ttest[icell, istim] = p_ttest_i

                evoked_cell = (stim_cell - base_cell) / (base_cell + 1e-7)
                evoked_i = np.mean(evoked_cell, axis=0) # trial averaged evoked resp
                evoked[icell, istim] = evoked_i
            
            # _, p_anova_cell = stats.f_oneway(np.array(base_cell_anova).flatten(), *stim_cell_anova)
            # p_anova[icell] = p_anova_cell
            # _, p_kruskal_cell = stats.kruskal(np.array(base_cell_anova).flatten(), *stim_cell_anova)
            # p_kruskal[icell] = p_kruskal_cell


    # cells responsive to image i: pass visually driven (anova OR kruskal) AND t-test AND amp threshold for *this* image
    p_sig = 0.05
    if area == 'LI':
        p_sig = 0.1 # relax sig for LI. used for vis_driven_ttest_bonferroni_jeff.mat
    evoked_thresh = 0.05
    img_driven = (p_ttest < p_sig) & (evoked > evoked_thresh)
    img_driven_msg = f'{img_driven.sum()} cells are image driven - with overlap between images, \n\
        proportion {np.round(img_driven.sum() / (ncell*nstim), 2)} out of {ncell*nstim} cell-stim combos. \n\
        1-N image evokes resp from {np.sum(img_driven, axis=0)} cells'
    logging.info(img_driven_msg)
    print(img_driven_msg)

    t = np.sum(img_driven, axis=1)
    logging.info(f'img driven cells are driven by {t[t>0]} images')

    # visually driven cells
    # NOTE changed to: driven by any image (pass t-test to any stim2, but with bonferroni correction) AND amp threshold. 
    loosen_sig = 1 # do not loosen, set ratio = 1
    p_bonferroni_corrected = loosen_sig * p_sig / nstim # number of possible stim2
    vis_driven = ((np.sum(p_ttest < p_bonferroni_corrected, axis=1) > 0) # any stim2 passes t-test
                & (np.sum(evoked > evoked_thresh, axis=1) > 0) # any stim2 evokes resp > thresh
                )
    vis_driven_msg = f'{vis_driven.sum()} cells are visually driven, \n\
        proportion {np.round(vis_driven.sum()/ncell, 2)} out of {ncell} cells'
    logging.info(vis_driven_msg)
    print(vis_driven_msg)

    # break

    ## save
    vis_driven_re = {}
    # vis_driven_re['p_anova'] = p_anova
    # vis_driven_re['p_kruskal'] = p_kruskal
    vis_driven_re['p_ttest'] = p_ttest
    vis_driven_re['evoked'] = evoked
    vis_driven_re['p_sig'] = p_sig
    vis_driven_re['p_bonferroni_corrected'] = p_bonferroni_corrected
    vis_driven_re['vis_driven'] = vis_driven
    vis_driven_re['img_driven'] = img_driven

    os.chdir(os.path.join(dir_inter, dir_sub))
    print(os.getcwd())
    savemat("vis_driven_ttest_bonferroni_jeff.mat", vis_driven_re)
    
    # with open('vis_driven_ttest_bonferroni_strict.pickle', 'wb') as f:
    #     pickle.dump(vis_driven_re, f)
    
    # break


# In[17]:


# session_mode = 'single session' # then irun must be appended to dir_sub !
# # session_mode = 'multi session' # then irun must be appended to dir_sub !

# ## set up logging
# plt.set_loglevel(level='warning') # turn off matplotlib debug logging
# logging.basicConfig(filename='C:/Users/ll357/Documents/inter/data/vis_driven.log', level=logging.DEBUG)
# logging.info('\n\n\n\n\n') # log several linebreaks to separate different runs
# logging.info(str(datetime.now()))

# for i in np.arange(nset): # TODO: this overwrites vis-driven.pickle multiple times for multisess data, fix later
#     # clear_output(wait=True)
#     logging.info(f'set {i+1} of {nset}')
#     mouse = meta.iloc[i].mouse.astype(str)
#     date = meta.iloc[i].date.astype(str)
#     area = meta.iloc[i].area
#     irun = '00' + meta.iloc[i].num.astype(int).astype(str)
#     logging.info(f'{mouse} {date} {irun} {area}')

#     ## check if resp_base_trialwise exists
#     ## load data
#     logging.info('waiting for mat')
#     print('waiting for mat')
#     while True:
#         try:
#             dir_sub = area + '_i' + mouse + '_' + date
#             if session_mode == 'single session':
#                 dir_sub += '_' + irun # else: multi sess, dont append irun
#             print(dir_sub)
#             try: # manual seg data has no suffix '_cellpose'
#                 dfof_trialwise = scipy.io.loadmat(os.path.join(dir_inter, dir_sub, 'resp_base_trialwise' + '.mat'))
#             except: # cellpose seg data has suffix '_cellpose'
#                 dir_sub += '_cellpose'
#                 dfof_trialwise = scipy.io.loadmat(os.path.join(dir_inter, dir_sub, 'resp_base_trialwise' + '.mat'))
#             logging.info('mat loaded')
#             break
#         except:
#             time.sleep(60)
#             clear_output(wait=True)
#             continue

#     if len(dfof_trialwise['dfof_ad_trial'].shape) == 2: # if contains only one ISI, likely 250
#         isi_mode = 'only one isi'
#         print(f'isi_mode is {isi_mode} (only one isi, grat_SF or bunny)')
#         print(dfof_trialwise['dfof_ad_trial'].shape, dfof_trialwise['dfof_ad_trial'][0,0].shape)
#         dfof_ad_trial = dfof_trialwise['dfof_ad_trial'] # do not subtract baseline here: stim resp will be compared to base resp
#         dfof_base_trial = dfof_trialwise['dfof_base_trial']

#     if len(dfof_trialwise['dfof_ad_trial'].shape) == 3: # if contains multiple ISI, likely inf-750-250
#         isi_mode = 'grat_8ori_3isi'
#         print(f'isi_mode is {isi_mode} (grating with 8 orientations and 3 ISI)')
#         print(dfof_trialwise['dfof_ad_trial'].shape, dfof_trialwise['dfof_ad_trial'][0,0,0].shape)
#         dfof_tg_trial = dfof_trialwise['dfof_tg_trial'][:,:,0] # only use the first ISI, aka no adapter trials
#         dfof_base2_trial = dfof_trialwise['dfof_base2_trial'][:,:,0]
#         print('to find visually driven cells and image driven cells, we should take only no-adapter trials and use target response. \n 8 possible target orientations will cover all needed cells')

#     ## stats
#     if isi_mode == 'only one isi':
#         print(isi_mode)

#         ncell = dfof_ad_trial.shape[0]
#         nstim = dfof_ad_trial.shape[1]

#         ## according to Ohki 2020
#         p_anova = np.ones((ncell, 1)) * np.nan # anova -> visually driven cells,  
#         p_kruskal = np.ones((ncell, 1)) * np.nan # t test -> certain image responsive cells,
#         p_ttest = np.ones((ncell, nstim)) * np.nan
#         evoked = np.ones((ncell, nstim)) * np.nan # amp thresh -> lower false positive rate

#         for icell in tqdm(np.arange(ncell)):
#             base_cell_anova = np.concatenate([dfof_base_trial[icell, stim] for stim in range(nstim)])
#             stim_cell_anova = [] 
#             for istim in np.arange(nstim):
#                 stim_cell_anova.append(np.array(dfof_ad_trial[icell, istim]).flatten())

#                 base_cell = dfof_base_trial[icell, istim]
#                 stim_cell = dfof_ad_trial[icell, istim]
#                 _, p_ttest_i = stats.ttest_ind(base_cell.flatten(), stim_cell.flatten(), equal_var=False, alternative='less')
#                 p_ttest[icell, istim] = p_ttest_i

#                 evoked_cell = (stim_cell - base_cell) / (base_cell + 1e-7)
#                 evoked_i = np.mean(evoked_cell, axis=0) # trial averaged evoked resp
#                 evoked[icell, istim] = evoked_i
            
#             _, p_anova_cell = stats.f_oneway(np.array(base_cell_anova).flatten(), *stim_cell_anova)
#             p_anova[icell] = p_anova_cell
#             _, p_kruskal_cell = stats.kruskal(np.array(base_cell_anova).flatten(), *stim_cell_anova)
#             p_kruskal[icell] = p_kruskal_cell

#     elif isi_mode == 'grat_8ori_3isi':
#         print(isi_mode)

#         ncell = dfof_tg_trial.shape[0]
#         nstim = dfof_tg_trial.shape[1]

#         ## according to Ohki 2020
#         p_anova = np.ones((ncell, 1)) * np.nan # anova -> visually driven cells,  
#         p_kruskal = np.ones((ncell, 1)) * np.nan # t test -> certain image responsive cells,
#         p_ttest = np.ones((ncell, nstim)) * np.nan
#         evoked = np.ones((ncell, nstim)) * np.nan # amp thresh -> lower false positive rate

#         for icell in tqdm(np.arange(ncell)):
#             base_cell_anova = np.concatenate([dfof_base2_trial[icell, stim] for stim in range(nstim)])
#             stim_cell_anova = [] 
#             for istim in np.arange(nstim):
#                 stim_cell_anova.append(np.array(dfof_tg_trial[icell, istim]).flatten())

#                 base_cell = dfof_base2_trial[icell, istim]
#                 stim_cell = dfof_tg_trial[icell, istim]
#                 _, p_ttest_i = stats.ttest_ind(base_cell.flatten(), stim_cell.flatten(), equal_var=False, alternative='less')
#                 p_ttest[icell, istim] = p_ttest_i

#                 evoked_cell = (stim_cell - base_cell) / (base_cell + 1e-7)
#                 evoked_i = np.mean(evoked_cell, axis=0) # trial averaged evoked resp
#                 evoked[icell, istim] = evoked_i
            
#             _, p_anova_cell = stats.f_oneway(np.array(base_cell_anova).flatten(), *stim_cell_anova)
#             p_anova[icell] = p_anova_cell
#             _, p_kruskal_cell = stats.kruskal(np.array(base_cell_anova).flatten(), *stim_cell_anova)
#             p_kruskal[icell] = p_kruskal_cell

#     # visually driven cells: pass (anova OR kruskal) AND amp threshold for >=1 image
#     p_sig = 0.05
#     vis_driven = ((p_anova < p_sig) | (p_kruskal < p_sig)) & (sum(evoked.T > 0.1) > 0).reshape(-1, 1)
#     logging.info(f'{vis_driven.sum()} cells are visually driven, \n\
#         proportion {np.round(vis_driven.sum()/ncell, 2)} out of {ncell} cells')

#     # cells responsive to image i: pass visually driven (anova OR kruskal) AND t-test AND amp threshold for *this* image
#     img_driven = vis_driven & (p_ttest < p_sig) & (evoked > 0.1)
#     logging.info(f'{img_driven.sum()} cells are image driven - with overlap between images, \n\
#         proportion {np.round(img_driven.sum() / (ncell*nstim), 2)} out of {ncell*nstim} cell-stim combos. \n\
#         1-N image evokes resp from {np.sum(img_driven, axis=0)} cells')

#     t = np.sum(img_driven, axis=1)
#     logging.info(f'img driven cells are driven by {t[t>0]} images')

#     ## save
#     vis_driven_re = {}
#     vis_driven_re['p_anova'] = p_anova
#     vis_driven_re['p_kruskal'] = p_kruskal
#     vis_driven_re['p_ttest'] = p_ttest
#     vis_driven_re['evoked'] = evoked
#     vis_driven_re['p_sig'] = p_sig
#     vis_driven_re['vis_driven'] = vis_driven
#     vis_driven_re['img_driven'] = img_driven

#     os.chdir(os.path.join(dir_inter, dir_sub))
#     print(os.getcwd())
#     with open('vis_driven.pickle', 'wb') as f:
#         pickle.dump(vis_driven_re, f)
    
#     # break


# # depre

# In[37]:


# # meta = pd.read_excel(mat_inter_path + 'adp_dataset_master.xlsx')

# # date_str = str(200721)
# # mouse_str = meta.loc[meta['date'] == int(date_str), 'mouse'].values#[0]
# # area_str = meta.loc[meta['date'] == int(date_str), 'area'].values[0]
# # if len(mouse_str) > 1:
# #     print('duplicate dates with maybe different mouse. select which mouse?')
# # else:
# #     mouse_str = str(mouse_str[0])
# # print(mouse_str, date_str)

# try:
#     dir_sub = area_str + '_i' + mouse_str + '_' + date_str + '_cellpose'
#     dfof_trialwise = scipy.io.loadmat(os.path.join(mat_inter_path, dir_sub, 'resp_base_trialwise' + '.mat')) # dir_sub[:-7] delete caiman in dir_sub
# except:
#     dir_sub = area_str + '_i' + mouse_str + '_' + date_str + ''
#     dfof_trialwise = scipy.io.loadmat(os.path.join(mat_inter_path, dir_sub, 'resp_base_trialwise' + '.mat')) # dir_sub[:-7] delete caiman in dir_sub
# if len(dfof_trialwise['dfof_ad_trial'].shape) == 2: # if contains only one ISI, likely 250
#     mode = 'only one isi'
#     print(f'mode is {mode} (only one isi, grat_SF or bunny)')
#     print(dfof_trialwise['dfof_ad_trial'].shape, dfof_trialwise['dfof_ad_trial'][0,0].shape)
#     dfof_ad_trial = dfof_trialwise['dfof_ad_trial'] # do not subtract baseline here: stim resp will be compared to base resp
#     dfof_base_trial = dfof_trialwise['dfof_base_trial']

# if len(dfof_trialwise['dfof_ad_trial'].shape) == 3: # if contains multiple ISI, likely inf-750-250
#     mode = 'grat_8ori_3isi'
#     print(f'mode is {mode} (grating with 8 orientations and 3 ISI)')
#     print(dfof_trialwise['dfof_ad_trial'].shape, dfof_trialwise['dfof_ad_trial'][0,0,0].shape)
#     dfof_tg_trial = dfof_trialwise['dfof_tg_trial'][:,:,0] # only use the first ISI, aka no adapter trials
#     dfof_base2_trial = dfof_trialwise['dfof_base2_trial'][:,:,0]
#     print('to find visually driven cells and image driven cells, we should take only no-adapter trials and use target response. \n 8 possible target orientations will cover all needed cells')


# In[38]:


# if mode == 'only one isi':
#     print(mode)

#     ncell = dfof_ad_trial.shape[0]
#     nstim = dfof_ad_trial.shape[1]

#     ## according to Ohki 2020
#     p_anova = np.ones((ncell, 1)) * np.nan # anova -> visually driven cells,  
#     p_kruskal = np.ones((ncell, 1)) * np.nan # t test -> certain image responsive cells,
#     p_ttest = np.ones((ncell, nstim)) * np.nan
#     evoked = np.ones((ncell, nstim)) * np.nan # amp thresh -> lower false positive rate

#     for icell in tqdm(np.arange(ncell)):
#         base_cell_anova = np.concatenate([dfof_base_trial[icell, stim] for stim in range(nstim)])
#         stim_cell_anova = [] 
#         for istim in np.arange(nstim):
#             stim_cell_anova.append(np.array(dfof_ad_trial[icell, istim]).flatten())

#             base_cell = dfof_base_trial[icell, istim]
#             stim_cell = dfof_ad_trial[icell, istim]
#             _, p_ttest_i = stats.ttest_ind(base_cell.flatten(), stim_cell.flatten(), equal_var=False, alternative='less')
#             p_ttest[icell, istim] = p_ttest_i

#             evoked_cell = (stim_cell - base_cell) / (base_cell + 1e-7)
#             evoked_i = np.mean(evoked_cell, axis=0) # trial averaged evoked resp
#             evoked[icell, istim] = evoked_i
        
#         _, p_anova_cell = stats.f_oneway(np.array(base_cell_anova).flatten(), *stim_cell_anova)
#         p_anova[icell] = p_anova_cell
#         _, p_kruskal_cell = stats.kruskal(np.array(base_cell_anova).flatten(), *stim_cell_anova)
#         p_kruskal[icell] = p_kruskal_cell

# elif mode == 'grat_8ori_3isi':
#     print(mode)

#     ncell = dfof_tg_trial.shape[0]
#     nstim = dfof_tg_trial.shape[1]

#     ## according to Ohki 2020
#     p_anova = np.ones((ncell, 1)) * np.nan # anova -> visually driven cells,  
#     p_kruskal = np.ones((ncell, 1)) * np.nan # t test -> certain image responsive cells,
#     p_ttest = np.ones((ncell, nstim)) * np.nan
#     evoked = np.ones((ncell, nstim)) * np.nan # amp thresh -> lower false positive rate

#     for icell in tqdm(np.arange(ncell)):
#         base_cell_anova = np.concatenate([dfof_base2_trial[icell, stim] for stim in range(nstim)])
#         stim_cell_anova = [] 
#         for istim in np.arange(nstim):
#             stim_cell_anova.append(np.array(dfof_tg_trial[icell, istim]).flatten())

#             base_cell = dfof_base2_trial[icell, istim]
#             stim_cell = dfof_tg_trial[icell, istim]
#             _, p_ttest_i = stats.ttest_ind(base_cell.flatten(), stim_cell.flatten(), equal_var=False, alternative='less')
#             p_ttest[icell, istim] = p_ttest_i

#             evoked_cell = (stim_cell - base_cell) / (base_cell + 1e-7)
#             evoked_i = np.mean(evoked_cell, axis=0) # trial averaged evoked resp
#             evoked[icell, istim] = evoked_i
        
#         _, p_anova_cell = stats.f_oneway(np.array(base_cell_anova).flatten(), *stim_cell_anova)
#         p_anova[icell] = p_anova_cell
#         _, p_kruskal_cell = stats.kruskal(np.array(base_cell_anova).flatten(), *stim_cell_anova)
#         p_kruskal[icell] = p_kruskal_cell


# In[45]:


# # visually driven cells: pass (anova OR kruskal) AND amp threshold for >=1 image
# p_sig = 0.05
# vis_driven = ((p_anova < p_sig) | (p_kruskal < p_sig)) & (sum(evoked.T > 0.1) > 0).reshape(-1, 1)
# print(f'{vis_driven.sum()} cells are visually driven, \n\
#     proportion {np.round(vis_driven.sum()/ncell, 2)} out of {ncell} cells')

# # cells responsive to image i: pass visually driven (anova OR kruskal) AND t-test AND amp threshold for *this* image
# img_driven = vis_driven & (p_ttest < p_sig) & (evoked > 0.1)
# print(f'{img_driven.sum()} cells are image driven - with overlap between images, \n\
#     proportion {np.round(img_driven.sum() / (ncell*nstim), 2)} out of {ncell*nstim} cell-stim combos. \n\
#     1-N image evokes resp from {np.sum(img_driven, axis=0)} cells')

# t = np.sum(img_driven, axis=1)
# print(f'img driven cells are driven by {t[t>0]} images')


# In[46]:


# vis_driven_re = {}
# vis_driven_re['p_anova'] = p_anova
# vis_driven_re['p_kruskal'] = p_kruskal
# vis_driven_re['p_ttest'] = p_ttest
# vis_driven_re['evoked'] = evoked
# vis_driven_re['p_sig'] = p_sig
# vis_driven_re['vis_driven'] = vis_driven
# vis_driven_re['img_driven'] = img_driven

# os.chdir(os.path.join(mat_inter_path, dir_sub))
# print(os.getcwd())

# with open('vis_driven.pickle', 'wb') as f:
#     pickle.dump(vis_driven_re, f)

