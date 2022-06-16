#!/usr/bin/env python
# coding: utf-8

# ## Prep

# In[44]:


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
meta = meta[meta.date == 210922].reset_index()
meta = meta[['mouse','date','area']]

meta.mouse = meta.mouse.astype(int)
meta.date = meta.date.astype(int).astype(str) + '_caiman'
meta = meta.head(1) # multisession, only keep one metadata
meta


# In[3]:


nset = len(meta.index); ncell = []; nori = 500; nisi = 1; nframe_trial = 137
dir_name = root_path + 'mat/'

vis_ad = np.empty([0,nori]); 
dfof_ad = np.empty([0,nori]); dfof_tg = np.empty([0,nori])
dfof_ad_std = np.empty([0,nori]); dfof_tg_std = np.empty([0,nori])
trace = np.empty([0,nori,nframe_trial])

for iset in np.arange(nset):
    dir_sub = str(meta.area[iset]) + '_i' + str(meta.mouse[iset]) + '_' + str(meta.date[iset])

    cell_prop = scipy.io.loadmat(os.path.join(dir_name, dir_sub, 'vis_driven' + '.mat'))
    dfof = scipy.io.loadmat(os.path.join(dir_name, dir_sub, 'dfof' + '.mat'))
    trace_align = scipy.io.loadmat(os.path.join(dir_name, dir_sub, 'trace_aligned' + '.mat'))
    
    ncell.append(len(cell_prop['vis_cell_ad']))    
    vis_ad = np.concatenate((vis_ad, cell_prop['sig_vis_ad']), axis=0)
    
    dfof_ad = np.concatenate((dfof_ad, dfof['dfof_ad']), axis=0)
    dfof_ad_std = np.concatenate((dfof_ad_std, dfof['dfof_ad_std']), axis=0)
    dfof_tg = np.concatenate((dfof_tg, dfof['dfof_tg']), axis=0)
    dfof_tg_std = np.concatenate((dfof_tg_std, dfof['dfof_tg_std']), axis=0)
    
    trace_flat = np.empty([ncell[iset],nori,nframe_trial]);
    for icell in np.arange(ncell[iset]):
        for iori in np.arange(nori):
            trace_flat[icell][iori][:] = trace_align['trace_avg'][icell][iori].flatten()
    trace = np.vstack((trace,trace_flat))

ncell, vis_ad.shape, dfof_ad.shape, dfof_tg.shape, trace.shape


# In[4]:


meta['ncell'] = ncell
meta


# In[5]:


mouse_cell = [item for item, count in zip(meta.mouse, meta.ncell) for i in range(count)]
area_cell = [item for item, count in zip(meta.area, meta.ncell) for i in range(count)]
meta_cell = pd.DataFrame({'mouse': mouse_cell, 'area': area_cell})
meta_cell.tail()


# ## Adaptation magnitude (todo: filter by vis_ad)

# adaptation mag = response to target ori==0 with adapter / response to adapter - 1  
# cell selection: vis_ad only, no dfof_ad thresholding

# In[6]:


adp_mag = dfof_tg / dfof_ad - 1
# adp_mag = adp_mag * vis_ad
adp_mag[adp_mag == 0] = np.nan
adp_mag.shape


# In[7]:


# plt.figure(figsize=(8,6))
# for iori in np.arange(nori):
#     plt.hist(adp_mag[:,iori], bins=np.linspace(-25,10,36), alpha=0.3)


# In[8]:


mag = adp_mag.flatten('F')
ad = dfof_ad.flatten('F')
tg = dfof_tg.flatten('F')
stim = [np.arange(nori)] * adp_mag.shape[0]
stim_flat = np.sort([item for sublist in stim for item in sublist])
stim_flat.shape, mag.shape, ad.shape


# In[9]:


df = pd.DataFrame({'stim':stim_flat, 'ad':ad, 'tg':tg, 'mag':mag, 'abs_mag':np.abs(mag)})
df


# In[10]:


mag_mean = df[['mag','stim']].groupby(['stim'], sort=False).mean().to_numpy().flatten()
mag_median = df[['mag','stim']].groupby(['stim'], sort=False).median().to_numpy().flatten()
mag_std = df[['mag','stim']].groupby(['stim'], sort=False).std().to_numpy().flatten()
mag_sem = df[['mag','stim']].groupby(['stim'], sort=False).sem().to_numpy().flatten()
df[['mag','stim']].groupby(['stim'], sort=False).describe()


# In[11]:


# plt.figure(figsize=(9,9))

# plt.errorbar(np.arange(nori), mag_median, yerr=mag_sem, zorder=100, 
#              color='silver', linewidth=3, capsize=3, capthick=3);
# sns.stripplot(x="stim", y="mag", data=df);

# # plt.ylim(-5,5);
# plt.ylabel('adaptation magnitude');
# plt.xlabel('stimulus identity');


# In[12]:


# plt.figure(figsize=(9,9))

# plt.errorbar(np.arange(nori), mag_mean, yerr=mag_sem, zorder=100, 
#              color='silver', linewidth=3, capsize=3, capthick=3);
# sns.stripplot(x="stim", y="mag", data=df);

# # plt.ylim(-5,5);
# plt.ylabel('adaptation magnitude');
# plt.xlabel('stimulus identity');


# ### adp_mag thresholding (todo: dfof_ad thresholding)

# In[13]:


plt.figure(figsize=(9,6))
ax = sns.scatterplot(data=df, x="ad", y="mag", hue="stim")
plt.xlim([-0.001,0.001]);
ax.set(title = 'adp mag changes with dfof_ad');


# cell selection: vis_ad only, with dfof_ad thresholding

# In[14]:


dfof_threshold = 0.0005
adp_threshold = 20 # 1, 2, 20

df_th = df.copy()
# df_th.loc[df_th[np.abs(df.ad) < dfof_threshold].index.to_numpy(),'mag'] = np.nan
df_th.loc[df_th[(df.ad) < dfof_threshold].index.to_numpy(),'mag'] = np.nan # not threshold by abs, bc 1-tail ttest originally?
df_th.loc[df_th[np.abs(df.mag) > adp_threshold].index.to_numpy(),'mag'] = np.nan
# bug / todo: fix vis_ad and filter normally by vis_ad, not filter by adp_mag

sns.histplot(data=df_th, x='mag');
print(df_th[np.isnan(df_th.mag)].shape[0] / df_th.shape[0])


# In[15]:


df_th[['mag','ad','tg','stim']].groupby(['stim'], sort=False).describe()


# In[16]:


mag_mean = df_th[['mag','stim']].groupby(['stim'], sort=False).mean().to_numpy().flatten()
mag_median = df_th[['mag','stim']].groupby(['stim'], sort=False).median().to_numpy().flatten()
mag_std = df_th[['mag','stim']].groupby(['stim'], sort=False).std().to_numpy().flatten()
mag_sem = df_th[['mag','stim']].groupby(['stim'], sort=False).sem().to_numpy().flatten()

# plt.figure(figsize=(9,9))
# plt.errorbar(np.arange(nori), mag_median, yerr=mag_sem, zorder=100, 
#              color='silver', linewidth=3, capsize=3, capthick=3);
# sns.stripplot(x="stim", y="mag", data=df_th);
# plt.ylabel('adaptation magnitude');
# plt.xlabel('stimulus identity');


# In[17]:


df_th[['mag','stim']].groupby(['stim'], sort=False).describe()


# In[18]:


mag_sort_id = mag_median.argsort()
mag_median_sorted = mag_median[mag_sort_id[::-1]]
mag_sem_sorted_med = mag_sem[mag_sort_id[::-1]]

mag_sort_id = mag_mean.argsort()
mag_mean_sorted = mag_mean[mag_sort_id[::-1]]
mag_sem_sorted_mean = mag_sem[mag_sort_id[::-1]]

mag_median_sorted[0], mag_median_sorted[-1], mag_mean_sorted[0], mag_mean_sorted[-1]


# ### real data, sorted adp

# In[19]:


plt.figure(figsize=(15,9))
plt.errorbar(np.arange(nori), mag_median_sorted, yerr=mag_sem_sorted_med, zorder=2, 
             color='silver', linewidth=3, capsize=3, capthick=3, alpha=0.5);
plt.errorbar(np.arange(nori), mag_mean_sorted, yerr=mag_sem_sorted_mean, zorder=1, 
             color='blue', linewidth=3, capsize=3, capthick=3, alpha=0.8);
plt.legend(['median', 'mean'])
plt.ylabel('adaptation magnitude');
plt.xlabel('stimulus identity sorted by adp mag');
plt.ylim(-2.3,0.1);
plt.title('real data');


# In[20]:


plt.figure(figsize=(15,9))
plt.plot(np.arange(nori), mag_median_sorted, zorder=2, 
             color='silver', linewidth=3, alpha=0.8);
plt.plot(np.arange(nori), mag_mean_sorted, zorder=1, 
             color='blue', linewidth=3, alpha=0.8);
plt.legend(['median', 'mean'])
plt.ylabel('adaptation magnitude');
plt.xlabel('stimulus identity sorted by adp mag');
plt.ylim(-2,0);


# ### shuffle
# ncell x nstim matrix shuffle, so that neuron population resp is not to a single stim, but random stims 

# In[21]:


adp_shuffle = adp_mag.copy()
print(adp_mag.shape) # ncell x nstim. should shuffle along rows (keep cell order, permute stim order)

adp_shuffle[np.abs(adp_shuffle) > adp_threshold] = np.nan

[np.random.shuffle(x) for x in adp_shuffle]; # shuffle every row


# In[22]:


mag = adp_shuffle.flatten('F')

ad = dfof_ad.flatten('F')
mag[np.abs(ad) < dfof_threshold] = np.nan

stim = [np.arange(nori)] * adp_mag.shape[0]
stim_flat = np.sort([item for sublist in stim for item in sublist])

df_shuffle = pd.DataFrame({'stim':stim_flat, 'mag':mag, 'abs_mag':np.abs(mag)})
df_shuffle.tail(5)


# In[23]:


mag_shuffle_mean = df_shuffle[['mag','stim']].groupby(['stim'], sort=False).mean().to_numpy().flatten()
mag_shuffle_median = df_shuffle[['mag','stim']].groupby(['stim'], sort=False).median().to_numpy().flatten()
mag_shuffle_std = df_shuffle[['mag','stim']].groupby(['stim'], sort=False).std().to_numpy().flatten()
mag_shuffle_sem = df_shuffle[['mag','stim']].groupby(['stim'], sort=False).sem().to_numpy().flatten()
df_shuffle[['mag','stim']].groupby(['stim'], sort=False).describe()


# In[24]:


mag_shuffle_sort_id = mag_shuffle_median.argsort()
mag_shuffle_median_sorted = mag_shuffle_median[mag_shuffle_sort_id[::-1]]
mag_shuffle_sem_sorted_med = mag_shuffle_sem[mag_shuffle_sort_id[::-1]]

mag_shuffle_sort_id = mag_shuffle_mean.argsort()
mag_shuffle_mean_sorted = mag_shuffle_mean[mag_shuffle_sort_id[::-1]]
mag_shuffle_sem_sorted_mean = mag_shuffle_sem[mag_shuffle_sort_id[::-1]]

mag_shuffle_median_sorted[0], mag_shuffle_median_sorted[-1], mag_shuffle_mean_sorted[0], mag_shuffle_mean_sorted[-1]


# In[25]:


plt.figure(figsize=(15,9))
plt.errorbar(np.arange(nori), mag_shuffle_median_sorted, yerr=mag_shuffle_sem_sorted_med, zorder=2, 
             color='silver', linewidth=3, capsize=3, capthick=3, alpha=0.5);
plt.errorbar(np.arange(nori), mag_shuffle_mean_sorted, yerr=mag_shuffle_sem_sorted_mean, zorder=1, 
             color='blue', linewidth=3, capsize=3, capthick=3, alpha=0.8);
plt.legend(['median', 'mean'])
plt.ylabel('adaptation mag_shufflenitude');
plt.xlabel('fake stimulus identity sorted by adp mag');
# plt.ylim(-2.3,0.1);
plt.title('shuffled stim');


# ### compare real vs shuffle

# In[26]:


plt.figure(figsize=(15,9))
plt.errorbar(np.arange(nori), mag_shuffle_median_sorted, yerr=mag_shuffle_sem_sorted_med, zorder=2, 
             color='red', linewidth=3, capsize=3, capthick=3, alpha=0.2);
plt.errorbar(np.arange(nori), mag_median_sorted, yerr=mag_sem_sorted_med, zorder=1, 
             color='blue', linewidth=3, capsize=3, capthick=3, alpha=0.2);
plt.legend(['shuffle', 'real'])
plt.ylabel('adaptation mag median');
plt.xlabel('stimulus identity sorted by adp mag');
plt.gca().set_ylim(top=0)


# In[27]:


plt.figure(figsize=(15,9))
plt.errorbar(np.arange(nori), mag_shuffle_mean_sorted, yerr=mag_shuffle_sem_sorted_mean, zorder=2, 
             color='red', linewidth=3, capsize=3, capthick=3, alpha=0.2);
plt.errorbar(np.arange(nori), mag_mean_sorted, yerr=mag_sem_sorted_mean, zorder=1, 
             color='blue', linewidth=3, capsize=3, capthick=3, alpha=0.2);
plt.legend(['shuffle', 'real'])
plt.ylabel('adaptation mag mean');
plt.xlabel('stimulus identity sorted by adp mag');
plt.gca().set_ylim(top=0)


# ## Visliz natural img sorted by adp

# In[28]:


# compare adp median vs mean distribution among stims

data = {"median": df_th[['mag','stim']].groupby(['stim'], sort=False).median().squeeze(), 
        "mean": df_th[['mag','stim']].groupby(['stim'], sort=False).mean().squeeze()}
ax = sns.histplot(data, bins=60, kde=True)

plt.axvline(df_th[['mag','stim']].groupby(['stim'], sort=False).median().mean()[0], color='b', alpha=0.7)
plt.axvline(df_th[['mag','stim']].groupby(['stim'], sort=False).mean().mean()[0], color='orange');


# In[113]:


# compare stim order sorted by median vs mean

stim_order_mean = df_th[['mag','stim']].groupby(['stim'], sort=False).mean().sort_values('mag', ascending=True).index.values
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


# In[86]:


# print img sorted by adp mean: cannot see any obvious feature trend

seq_len = 500
stim_top = stim_order_mean[:seq_len]
stim_bottom = stim_order_mean[-seq_len:]
stim_folder = root_path + 'code/mwork/cut/keep_500/'
save_path = root_path + 'code/mwork/cut/'

for idx, stim in enumerate(stim_top):
    if stim == 0:
        continue
    img_name = stim_folder + str(stim) + '.jpg'
    img = Image.open(img_name)
#     print(stim)
#     display(img)
    adp_stim = adp_mean[idx]    
    img.save(save_path + f'adp {adp_stim:.2f} stim {stim}.png')


# In[114]:


seq_len = 10
stim_top = stim_order_mean[:seq_len]
stim_bottom = stim_order_mean[-seq_len:]
stim_mid = stim_order_mean[nori//2-seq_len//2 : nori//2+seq_len//2]

# np.argwhere(stim_order_mean == stim_top[0]), np.argwhere(stim_order_mean == stim_top[-1])
np.argwhere(stim_order_mean == stim_mid[0]), np.argwhere(stim_order_mean == stim_mid[-1])
# np.argwhere(stim_order_mean == stim_bottom[0]), np.argwhere(stim_order_mean == stim_bottom[-1])


# ### construct bunnytop stim x30

# In[116]:


stim_top, stim_mid, stim_bottom


# In[119]:


plt.plot(np.sort(stim_top))
plt.plot(np.sort(stim_mid))
plt.plot(np.sort(stim_bottom))


# In[136]:


# copy & rename top, mid, bottom imgs to 1-30.jpg to satisfy mwork

img_old = [462,318,408,179,69,13,28,8,301,65,    435,422,429,81,267,402,1,74,136,19,    494,47,303,187,175,319,413,271,210,115]
img_new = list(np.arange(31)[1:])

for i in np.arange(len(img_old)):
    img_id = img_old[i]
    img_id2 = img_new[i]
    src = root_path + f'code/mwork/bunny500/Image/{img_id}.jpg'
    dst = root_path + f'code/mwork/bunnytop/Image/{img_id2}.jpg'
    shutil.copyfile(src, dst, follow_symlinks=True)


# ## Trace (cannot see bc no vis_ad filter)

# In[87]:


trace.shape, trace_mean.shape


# In[74]:


trace_mean = np.mean(trace, axis=0)
trace_std = np.std(trace, axis=0)
trace_sem = trace_std / np.sqrt(trace.shape[0])


# In[85]:


# fig, ax = plt.subplots(figsize=(20,15))
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
plt.ylim(-0.02,0.10)
plt.xlabel('frame number')
plt.ylabel('dfof')
plt.show()


# In[ ]:




