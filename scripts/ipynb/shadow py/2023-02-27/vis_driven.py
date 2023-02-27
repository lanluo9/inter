#!/usr/bin/env python
# coding: utf-8

# In[35]:


import numpy as np
import pandas as pd
import scipy.io
from scipy import stats

# import seaborn as sns
import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

from tqdm import tqdm
import os
import pickle


# In[36]:


mat_inter_path = r'\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_Staff\home\lan\Data\2P_images\mat_inter/'.replace('\\', '/')
meta = pd.read_excel(mat_inter_path + 'adp_dataset_master.xlsx')

date_str = str(220529)
mouse_str = meta.loc[meta['date'] == int(date_str), 'mouse'].values#[0]
if len(mouse_str) > 1:
    print('duplicate dates with maybe different mouse. select which mouse?')
else:
    mouse_str = str(mouse_str[0])
print(mouse_str, date_str)

dir_sub = 'V1' + '_i' + mouse_str + '_' + date_str + '_lindsey' # '_cellpose'
dfof_trialwise = scipy.io.loadmat(os.path.join(mat_inter_path, dir_sub, 'resp_base_trialwise' + '.mat')) # dir_sub[:-7] delete caiman in dir_sub
print(dfof_trialwise['dfof_ad_trial'].shape, dfof_trialwise['dfof_ad_trial'][0,0].shape)

dfof_ad_trial = dfof_trialwise['dfof_ad_trial'] # do not subtract baseline here: stim resp will be compared to base resp
dfof_tg_trial = dfof_trialwise['dfof_tg_trial']


# # AND gate for cell filter
# according to Ohki 2020  
# anova -> visually driven cells,  
# t test -> certain image responsive cells,  
# amp thresh -> lower false positive rate

# In[37]:


ncell = dfof_ad_trial.shape[0]
nstim = dfof_ad_trial.shape[1]

p_anova = np.ones((ncell, 1)) * np.nan
p_kruskal = np.ones((ncell, 1)) * np.nan
p_ttest = np.ones((ncell, nstim)) * np.nan
evoked = np.ones((ncell, nstim)) * np.nan

for icell in tqdm(np.arange(ncell)):
    base_cell_anova = np.concatenate([dfof_trialwise['dfof_base_trial'][icell, stim] for stim in range(nstim)])
    stim_cell_anova = [] 
    for istim in np.arange(nstim):
        stim_cell_anova.append(np.array(dfof_ad_trial[icell, istim]).flatten())

        base_cell = dfof_trialwise['dfof_base_trial'][icell, istim]
        stim_cell = dfof_trialwise['dfof_ad_trial'][icell, istim]
        _, p_ttest_i = stats.ttest_ind(base_cell.flatten(), stim_cell.flatten(), equal_var=False, alternative='less')
        p_ttest[icell, istim] = p_ttest_i

        evoked_cell = (stim_cell - base_cell) / (base_cell + 1e-7)
        evoked_i = np.mean(evoked_cell, axis=0) # trial averaged evoked resp
        evoked[icell, istim] = evoked_i
    
    _, p_anova_cell = stats.f_oneway(np.array(base_cell_anova).flatten(), *stim_cell_anova)
    p_anova[icell] = p_anova_cell
    _, p_kruskal_cell = stats.kruskal(np.array(base_cell_anova).flatten(), *stim_cell_anova)
    p_kruskal[icell] = p_kruskal_cell


# In[38]:


# visually driven cells: pass (anova OR kruskal) AND amp threshold for >=1 image
p_sig = 0.01
vis_driven = ((p_anova < p_sig) | (p_kruskal < p_sig)) & (sum(evoked.T > 0.1) > 0).reshape(-1, 1)
print(f'{vis_driven.sum()} cells are visually driven, \n    proportion {np.round(vis_driven.sum()/ncell, 2)} out of {ncell} cells')

# cells responsive to image i: pass visually driven (anova OR kruskal) AND t-test AND amp threshold for *this* image
img_driven = vis_driven & (p_ttest < p_sig) & (evoked > 0.1)
print(f'{img_driven.sum()} cells are image driven - with overlap between images, \n    proportion {np.round(img_driven.sum() / (ncell*nstim), 2)} out of {ncell*nstim} cell-stim combos. \n    1-30 image evokes resp from {np.sum(img_driven, axis=0)} cells')

t = np.sum(img_driven, axis=1)
print(f'img driven cells are driven by {t[t>0]} images')


# In[39]:


os.chdir(os.path.join(mat_inter_path, dir_sub))
print(os.getcwd())

# with open('vis_driven' + '.pickle', 'wb') as handle:
#     pickle.dump(p_anova, handle, protocol=pickle.HIGHEST_PROTOCOL)
#     pickle.dump(p_kruskal, handle, protocol=pickle.HIGHEST_PROTOCOL)
#     pickle.dump(evoked, handle, protocol=pickle.HIGHEST_PROTOCOL)
#     pickle.dump(p_ttest, handle, protocol=pickle.HIGHEST_PROTOCOL)

# with open('vis_driven.pickle', 'wb') as f:
#     pickle.dump([p_anova, p_kruskal, evoked, p_ttest], f)

vis_driven_re = {}
vis_driven_re['p_anova'] = p_anova
vis_driven_re['p_kruskal'] = p_kruskal
vis_driven_re['evoked'] = evoked
vis_driven_re['p_ttest'] = p_ttest
vis_driven_re['vis_driven'] = vis_driven
vis_driven_re['img_driven'] = img_driven

with open('vis_driven.pickle', 'wb') as f:
    pickle.dump(vis_driven_re, f)

# with open('vis_driven.pickle', 'rb') as f: 
#     p_anova, p_kruskal, evoked, p_ttest = pickle.load(f)

