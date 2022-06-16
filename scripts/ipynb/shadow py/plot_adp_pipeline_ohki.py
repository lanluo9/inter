#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import pandas as pd
import scipy.io
from sklearn.linear_model import RidgeCV

import seaborn as sns
import matplotlib.pyplot as plt
from ipywidgets import interactive
get_ipython().run_line_magic('matplotlib', 'inline')


# ## Load data

# In[2]:


meta = pd.read_excel('C:/Users/lan/Documents/repos/inter/mat/adp_dataset_master.xlsx', index_col=None)
meta = meta[meta.seg == 'segmented']
meta = meta[['mouse','date','area']]


# In[11]:


nset = len(meta.index); ncell = []; nori = 8; nisi = 3; 
dir_name = 'C:\\Users\\lan\\Documents\\repos\\inter\\mat\\'

vis_driven = np.empty([0,1]); ori_driven = np.empty([0,nori]); ori_driven_cell = np.empty([0,1]);
dfof_tg = np.empty([0,nori,nisi]); dfof_tg_std = dfof_tg; dfof_tg_sem = dfof_tg

for iset in np.arange(nset):
    dir_sub = str(meta.area[iset]) + '_i' + str(meta.mouse[iset]) + '_' + str(meta.date[iset]) + '_ohki'

    cell_prop = scipy.io.loadmat(os.path.join(dir_name, dir_sub, 'cell_property.mat'))
    dfof = scipy.io.loadmat(os.path.join(dir_name, dir_sub, 'dfof_tg' + '.mat'))
    
    ncell.append(len(cell_prop['vis_driven']))    
    vis_driven = np.concatenate((vis_driven, cell_prop['vis_driven']), axis=0)
    ori_driven = np.concatenate((ori_driven, cell_prop['ori_driven']), axis=0)
    ori_driven_cell = np.concatenate((ori_driven_cell, cell_prop['ori_driven_cell']), axis=0)
    
    dfof_tg = np.concatenate((dfof_tg, dfof['dfof_tg']), axis=0)
    dfof_tg_std = np.concatenate((dfof_tg_std, dfof['dfof_tg_std']), axis=0)
    dfof_tg_sem = np.concatenate((dfof_tg_sem, dfof['dfof_tg_sem']), axis=0)

ncell, vis_driven.shape, ori_driven.shape, ori_driven_cell.shape, dfof_tg.shape, dfof_tg_std.shape, dfof_tg_sem.shape


# In[13]:


meta['ncell'] = ncell
# meta = meta.replace({'area' : { 'V1':1, 'LM':2, 'LI':3 }})
mouse_cell = [item for item, count in zip(meta.mouse, meta.ncell) for i in range(count)]
area_cell = [item for item, count in zip(meta.area, meta.ncell) for i in range(count)]
meta_cell = pd.DataFrame({'mouse': mouse_cell, 'area': area_cell})
meta_cell 


# ## Ori responsive distribution
# to determine cell selection in training set

# In[228]:


df = meta_cell.copy()
df['vis_driven'] = vis_driven.flatten()

attach = pd.DataFrame(columns = ['ori_driven'])
for row in ori_driven:
    temp = pd.DataFrame(columns = ['ori_driven']); temp['ori_driven'] = [row]
    attach = pd.concat([attach, temp])

df2 = pd.concat([df, attach.reset_index()], axis=1)
df2.drop(['index'], axis=1, inplace=True)

df3 = df2.copy()
df3 = df3[df3.vis_driven == 1].reset_index().drop(['index'], axis=1)
df3


# In[343]:


x = ['V1', 'LM', 'LI']
y = np.flip(np.asarray(df2[['area','vis_driven']].groupby('area').mean()).flatten())
plt.bar(x, y)
plt.ylim([0,1])
plt.show()


# In[302]:


ori_driven_area = np.empty([len(np.unique(df3.area)), nori]); iarea = 0
for area in np.unique(df3.area):
    df_area = df3[df3.area == area]
    df_area = df_area.reset_index()
    t = df_area['ori_driven'].to_numpy()
    t = np.concatenate(t).ravel().reshape([-1, nori])
    ori_driven_area[iarea,:] = np.sum(t, axis=0)
    iarea += 1
np.unique(df3.area), ori_driven_area


# In[348]:


ori_percent = ori_driven_area / ncell_area
ori_percent_area = np.flip(np.mean(ori_percent, axis=1, keepdims = 1))

x = ['V1', 'LM', 'LI']
y = ori_percent_area.flatten()
plt.bar(x, y)
plt.ylim([0,1])
plt.show()


# In[349]:


ori_percent_area


# In[340]:


plt.figure(figsize=(15, 10))
ax = plt.subplot(121, projection='polar')
N = nori * 2
theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
width = np.pi / nori
color_list = ['C0', 'C1', 'C3']
for i in np.flip(np.arange(len(np.unique(df3.area)))):
    y = np.concatenate((ori_driven_area[i], ori_driven_area[i]))
    ax.bar(theta, y, width=width, bottom=0.0, color=color_list[i], alpha=1-0.35*i)
plt.yticks(np.arange(0, 90, step=15))
ax.legend(['V1', 'LM', 'LI'])

ncell_area = np.asarray(df2[['area','ori_driven']].groupby(['area']).count())
ax = plt.subplot(122, projection='polar')
for i in np.flip(np.arange(len(np.unique(df3.area)))):
    y = np.concatenate((ori_driven_area[i]/ncell_area[i], ori_driven_area[i]/ncell_area[i]))
    y = 100*y # conver to %
    ax.bar(theta, y, width=width, bottom=0.0, color=color_list[i], alpha=0.2)
# plt.yticks(np.arange(0, 21, step=3))
ax.legend(['V1', 'LM', 'LI'])
plt.show()


# In[237]:


df.groupby(['area','mouse']).describe()


# ## Fit encoding model

# In[438]:


nset = len(meta.index);
dir_name = 'C:/Users/lan/Documents/repos/inter/mat/'
R2 = np.empty([0])

for iset in np.arange(nset):
    dir_sub = str(meta.area[iset]) + '_i' + str(meta.mouse[iset]) + '_' + str(meta.date[iset]) + '_ohki'
    dfof_tg_noad_z_score = scipy.io.loadmat(os.path.join(dir_name, dir_sub, 'dfof_tg_noad_z_score.mat'))
    feature_trial = scipy.io.loadmat(os.path.join(dir_name, dir_sub, 'feature_trial.mat'))
    
    z_score = dfof_tg_noad_z_score['z_score']
    F_trial = feature_trial['F_trial']
    ncell_set = z_score.shape[0]
    
    R2_set = np.zeros([ncell_set])
    for icell in np.arange(ncell_set):
        X = F_trial.T; y = z_score[icell,:]
        clf = RidgeCV(alphas=[1e-3, 1e-2, 1e-1, 1], cv=10, fit_intercept=True).fit(X, y)
        R2_set[icell] = clf.score(X, y)
    R2 = np.concatenate((R2, R2_set), axis=0)
    
R2.shape


# In[445]:


df2['R2'] = R2
df2['ori_driven_cell'] = ori_driven_cell
df4 = df2.copy()
df4 = df4[(df4.vis_driven == 1) & (df4.ori_driven_cell == 1)].reset_index().drop(['index'], axis=1)
df4


# In[443]:


def violinplot(df, ix, iy, ihue, iylabel, ititle, fig_size=[8,5]):
    # ivar = string
    # violinplot(df_adp_fano, "area", "adp_fano", "isi", 'fano factor change after adp', 'adaptation impacts variability')
    
    sns.set_style("whitegrid")
    plt.figure(figsize=(fig_size))
    ax = sns.violinplot(data=df, x=ix, y=iy, hue=ihue, split=True, inner="quart", palette="Set3")
    sns.despine(left=True)
    ax.set(ylabel = iylabel, title = ititle)
    plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# In[450]:


sns.set_style("whitegrid")
plt.figure(figsize=[8,5])
ax = sns.violinplot(data=df4, x="area", y="R2", inner="quart", palette="Set3")
sns.despine(left=True)
ax.set(ylabel = 'encoding model R2', title = '')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# In[451]:


df4[['area','R2']].groupby(['area'], sort=False).describe()


# In[385]:


# clf.coef_.shape, clf.intercept_, clf.alpha_, clf.best_score_, clf.score(X, y)

