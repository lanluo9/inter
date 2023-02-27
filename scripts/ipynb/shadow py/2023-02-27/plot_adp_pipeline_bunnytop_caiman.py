#!/usr/bin/env python
# coding: utf-8

# ## Prep

# In[1]:


import numpy as np
import pandas as pd
import scipy.io

import seaborn as sns
import matplotlib.pyplot as plt
# from ipywidgets import interactive
get_ipython().run_line_magic('matplotlib', 'inline')

from tqdm import tqdm
import os
import shutil
from PIL import Image


# In[2]:


root_path = 'C:/Users/ll357/Documents/inter/'
meta = pd.read_excel(root_path + 'mat/adp_dataset_master.xlsx', index_col=None)
meta = meta[meta.date == 211222].reset_index()
meta = meta[['mouse','date','area']]

meta.mouse = meta.mouse.astype(int)
meta.date = meta.date.astype(int).astype(str) + '_caiman'
meta = meta.head(1) # multisession, only keep one metadata
meta


# In[14]:


nset = len(meta.index); ncell = []; nori = 30; nisi = 1; nframe_trial = 143
dir_name = root_path + 'mat/'

# vis_ad = np.empty([0,nori]); 
dfof_ad = np.empty([0,nori]); dfof_tg = np.empty([0,nori])
dfof_ad_std = np.empty([0,nori]); dfof_tg_std = np.empty([0,nori])
trace = np.empty([0,nori,nframe_trial])

for iset in np.arange(nset):
    dir_sub = str(meta.area[iset]) + '_i' + str(meta.mouse[iset]) + '_' + str(meta.date[iset])

    dfof = scipy.io.loadmat(os.path.join(dir_name, dir_sub, 'dfof' + '.mat'))
    trace_align = scipy.io.loadmat(os.path.join(dir_name, dir_sub, 'trace_aligned' + '.mat'))
    cell_prop = scipy.io.loadmat(os.path.join(dir_name, dir_sub, 'vis_driven' + '.mat'))
    
    ncell.append(len(cell_prop['vis_cell_ad']))    
#     vis_ad = np.concatenate((vis_ad, cell_prop['sig_vis_ad']), axis=0)
    
    dfof_ad = np.concatenate((dfof_ad, dfof['dfof_ad']), axis=0)
    dfof_ad_std = np.concatenate((dfof_ad_std, dfof['dfof_ad_std']), axis=0)
    dfof_tg = np.concatenate((dfof_tg, dfof['dfof_tg']), axis=0)
    dfof_tg_std = np.concatenate((dfof_tg_std, dfof['dfof_tg_std']), axis=0)
    
    trace_flat = np.empty([ncell[iset],nori,nframe_trial]);
    for icell in np.arange(ncell[iset]):
        for iori in np.arange(nori):
            trace_flat[icell][iori][:] = trace_align['trace_avg'][icell][iori].flatten()
    trace = np.vstack((trace,trace_flat))

ncell, dfof_ad.shape, dfof_tg.shape, trace.shape #vis_ad.shape,


# In[15]:


meta['ncell'] = ncell
meta


# In[16]:


mouse_cell = [item for item, count in zip(meta.mouse, meta.ncell) for i in range(count)]
area_cell = [item for item, count in zip(meta.area, meta.ncell) for i in range(count)]
meta_cell = pd.DataFrame({'mouse': mouse_cell, 'area': area_cell})
meta_cell.tail()


# ## Adaptation magnitude (todo: filter by vis_ad)

# adaptation mag = response to target ori==0 with adapter / response to adapter - 1  
# cell selection: vis_ad only, no dfof_ad thresholding

# In[17]:


adp_mag = dfof_tg / dfof_ad - 1
# adp_mag = adp_mag * vis_ad
adp_mag[adp_mag == 0] = np.nan
adp_mag.shape


# In[19]:


# plt.figure(figsize=(8,6))
# for iori in np.arange(nori):
#     plt.hist(adp_mag[:,iori], bins=np.linspace(-25,10,36), alpha=0.3)


# In[18]:


mag = adp_mag.flatten('F')
ad = dfof_ad.flatten('F')
tg = dfof_tg.flatten('F')
stim = [np.arange(nori)] * adp_mag.shape[0]
stim_flat = np.sort([item for sublist in stim for item in sublist])
stim_flat.shape, mag.shape, ad.shape


# In[20]:


df = pd.DataFrame({'stim':stim_flat, 'ad':ad, 'tg':tg, 'mag':mag, 'abs_mag':np.abs(mag)})
# df


# In[21]:


mag_mean = df[['mag','stim']].groupby(['stim'], sort=False).mean().to_numpy().flatten()
mag_median = df[['mag','stim']].groupby(['stim'], sort=False).median().to_numpy().flatten()
mag_std = df[['mag','stim']].groupby(['stim'], sort=False).std().to_numpy().flatten()
mag_sem = df[['mag','stim']].groupby(['stim'], sort=False).sem().to_numpy().flatten()
# df[['mag','stim']].groupby(['stim'], sort=False).describe()


# ### adp_mag thresholding (todo: dfof_ad thresholding)

# In[26]:


plt.figure(figsize=(9,6))
ax = sns.scatterplot(data=df, x="ad", y="mag", hue="stim")
plt.xlim([-0.001,0.001]);
ax.set(title = 'adp mag changes with dfof_ad');


# cell selection: vis_ad only, with dfof_ad thresholding

# In[121]:


sns.histplot(data=df_th, x='ad');
sum(df_th.ad < 0) / len(df_th)


# In[55]:


dfof_threshold = 0.00025
adp_threshold = 20 # 1, 2, 20

df_th = df.copy()
# df_th.loc[df_th[np.abs(df.ad) < dfof_threshold].index.to_numpy(),'mag'] = np.nan
df_th.loc[df_th[(df.ad) < dfof_threshold].index.to_numpy(),'mag'] = np.nan # not threshold by abs, bc 1-tail ttest originally?
df_th.loc[df_th[np.abs(df.mag) > adp_threshold].index.to_numpy(),'mag'] = np.nan
# bug / todo: fix vis_ad and filter normally by vis_ad, not filter by adp_mag

sns.histplot(data=df_th, x='mag');
print(df_th[np.isnan(df_th.mag)].shape[0] / df_th.shape[0])


# ### groupby bunnytop 3 groups

# In[56]:


df_th['group'] = ''
df_th.loc[df_th.stim < 30, 'group'] = 'bottom'
df_th.loc[df_th.stim < 20, 'group'] = 'mid'
df_th.loc[df_th.stim < 10, 'group'] = 'top' # sorted by ascending adp, "top" adp imgs are the most adapted (lowest neg value)
df_th


# In[124]:


mag_mean = df_th[['mag','group']].groupby(['group'], sort=False).mean().to_numpy().flatten()
mag_median = df_th[['mag','group']].groupby(['group'], sort=False).median().to_numpy().flatten()
mag_std = df_th[['mag','group']].groupby(['group'], sort=False).std().to_numpy().flatten()
mag_sem = df_th[['mag','group']].groupby(['group'], sort=False).sem().to_numpy().flatten()

ngroup = len(np.unique(df_th.group))
plt.figure(figsize=(9,9))
plt.errorbar(np.arange(len(np.unique(df_th.group))), mag_median, yerr=mag_sem, zorder=100, 
             color='silver', linewidth=3, capsize=3, capthick=3);
plt.errorbar(np.arange(len(np.unique(df_th.group))), mag_mean, yerr=mag_sem, zorder=100, 
             color='blue', linewidth=3, capsize=3, capthick=3);
# sns.stripplot(x="group", y="mag", data=df_th);
plt.ylabel('adaptation magnitude');
plt.xlabel('stimulus identity');
plt.xticks([0,1,2], ['top','mid','bottom'])
# plt.ylim(-1.1,0)
plt.gca().set_ylim(top=0);


# In[65]:


df_th[['mag','group']].groupby(['group'], sort=False).describe()


# In[125]:


# mag_sort_id = mag_median.argsort()
# mag_median_sorted = mag_median[mag_sort_id[::-1]]
# mag_sem_sorted_med = mag_sem[mag_sort_id[::-1]]

# mag_sort_id = mag_mean.argsort()
# mag_mean_sorted = mag_mean[mag_sort_id[::-1]]
# mag_sem_sorted_mean = mag_sem[mag_sort_id[::-1]]

# mag_median_sorted[0], mag_median_sorted[-1], mag_mean_sorted[0], mag_mean_sorted[-1]


# ### one way anova btw 3 groups

# In[106]:


from scipy.stats import f_oneway
f_oneway(df_th.loc[(df_th.group == 'top') & (~np.isnan(df_th.mag))].mag.to_numpy(),
         df_th.loc[(df_th.group == 'mid') & (~np.isnan(df_th.mag))].mag.to_numpy(),
         df_th.loc[(df_th.group == 'bottom') & (~np.isnan(df_th.mag))].mag.to_numpy())


# ## compare img order sorted by adp
# ### bunny500 vs bunnytop

# In[80]:


# compare adp median vs mean distribution among stims
data = {"median": df_th[['mag','stim']].groupby(['stim'], sort=False).median().squeeze(), 
        "mean": df_th[['mag','stim']].groupby(['stim'], sort=False).mean().squeeze()}
ax = sns.histplot(data, bins=60, kde=True)

plt.axvline(df_th[['mag','stim']].groupby(['stim'], sort=False).median().mean()[0], color='b', alpha=0.7)
plt.axvline(df_th[['mag','stim']].groupby(['stim'], sort=False).mean().mean()[0], color='orange');


# In[81]:


# compare stim order sorted by median vs mean

stim_order_mean = df_th[['mag','stim']].groupby(['stim'], sort=False).mean().sort_values('mag').index.values
adp_mean = np.sort(mag_mean.copy())
stim_order_median = df_th[['mag','stim']].groupby(['stim'], sort=False).median().sort_values('mag').index.values
adp_median = np.sort(mag_median.copy())

seq_len = 10
set1 = set(stim_order_mean[:seq_len])
set2 = set(stim_order_median[:seq_len])
intersect_ratio = len(set1.intersection(set2)) / seq_len
print(intersect_ratio)

set1 = set(stim_order_mean[-seq_len:])
set2 = set(stim_order_median[-seq_len:])
intersect_ratio = len(set1.intersection(set2)) / seq_len
print(intersect_ratio)

intersect_ratio = 0
count = 0
for i in np.arange(10000):
    set1 = set(np.random.choice(nori, seq_len, replace=False)) # nori = nstim
    set2 = set(np.random.choice(nori, seq_len, replace=False))
    intersect_ratio = intersect_ratio + len(set1.intersection(set2)) / seq_len
    count = count + 1
intersect_ratio = intersect_ratio / count
print(intersect_ratio)


# In[82]:


stim_order_mean


# In[84]:


np.arange(31)[1:]


# In[85]:


plt.plot(stim_order_mean)
plt.plot(np.arange(31)[1:])
# plt.plot(np.sort(stim_bottom))


# ### bunnytop sess1 vs 3

# In[ ]:





# ## Trace (cannot see bc no vis_ad filter)

# In[87]:


trace.shape#, trace_mean.shape


# In[88]:


trace_mean = np.mean(trace, axis=0)
trace_std = np.std(trace, axis=0)
trace_sem = trace_std / np.sqrt(trace.shape[0])


# In[90]:


fig, ax = plt.subplots(figsize=(9,6))
for iori in np.arange(nori):
    ax.plot(np.arange(trace.shape[2]), trace_mean[iori,:].flatten(), )
#             color=sns.color_palette("bright")[iori])
    ax.fill_between(np.arange(trace.shape[2]), 
                    trace_mean[iori,:].flatten() + trace_sem[iori,:].flatten(), 
                    trace_mean[iori,:].flatten() - trace_sem[iori,:].flatten(),
                    alpha=0.1) # color=sns.color_palette("bright")[iori], 
plt.grid('minor')
plt.xlim(0, 50)
plt.xlabel('frame number')
plt.ylabel('dfof')
plt.show()


# ## Vis driven
# ### 1 way ANOVA

# In[118]:


baseline = trace[:,:,:4]
baseline = np.mean(np.mean(baseline, 2),1)
baseline.shape


# In[111]:


dfof_ad.shape

