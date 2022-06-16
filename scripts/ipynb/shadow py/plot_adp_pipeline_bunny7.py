#!/usr/bin/env python
# coding: utf-8

# In[3]:


import os
import numpy as np
import pandas as pd
import scipy.io

import seaborn as sns
import matplotlib.pyplot as plt
# from ipywidgets import interactive
get_ipython().run_line_magic('matplotlib', 'inline')

from tqdm import tqdm
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.linear_model import LinearRegression
import random


# ## Load data

# In[6]:


root_path = 'C:/Users/ll357/Documents/inter/'
meta = pd.read_excel(root_path + 'mat/adp_dataset_master.xlsx', index_col=None)
meta = meta[meta.date == 210616].reset_index()
meta = meta[['mouse','date','area']]

meta.mouse = meta.mouse.astype(int)
meta.date = meta.date.astype(int)
meta


# In[8]:


nset = len(meta.index); ncell = []; nori = 7; nisi = 1; nframe_trial = 77
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


# In[9]:


meta['ncell'] = ncell
meta


# In[10]:


mouse_cell = [item for item, count in zip(meta.mouse, meta.ncell) for i in range(count)]
area_cell = [item for item, count in zip(meta.area, meta.ncell) for i in range(count)]
meta_cell = pd.DataFrame({'mouse': mouse_cell, 'area': area_cell})
meta_cell.tail()


# ## Adaptation magnitude

# adaptation mag = response to target ori==0 with adapter / response to adapter - 1  
# cell selection: vis_ad only, no dfof_ad thresholding

# In[11]:


adp_mag = dfof_tg / dfof_ad - 1
adp_mag = adp_mag * vis_ad
adp_mag[adp_mag == 0] = np.nan
adp_mag.shape


# In[12]:


plt.figure(figsize=(8,6))
for iori in np.arange(nori):
#     plt.hist(adp_mag[:,iori], alpha=0.3)
    plt.hist(adp_mag[:,iori], bins=np.linspace(-25,10,36), alpha=0.3)


# In[13]:


plt.figure(figsize=(8,6))
for iori in np.arange(nori):
#     plt.hist(adp_mag[:,iori], alpha=0.3)
    plt.hist(adp_mag[:,iori], bins=np.linspace(-25,10,36), alpha=0.3)
plt.xlim(-2,2)


# In[14]:


mag = adp_mag.flatten('F')
ad = dfof_ad.flatten('F')
tg = dfof_tg.flatten('F')
stim = [np.arange(nori)] * adp_mag.shape[0]
stim_flat = np.sort([item for sublist in stim for item in sublist])
stim_flat.shape, mag.shape, ad.shape


# In[15]:


df = pd.DataFrame({'stim':stim_flat, 'ad':ad, 'tg':tg, 'mag':mag, 'abs_mag':np.abs(mag)})
df


# In[17]:


mag_mean = df[['mag','stim']].groupby(['stim'], sort=False).mean().to_numpy().flatten()
mag_median = df[['mag','stim']].groupby(['stim'], sort=False).median().to_numpy().flatten()
mag_std = df[['mag','stim']].groupby(['stim'], sort=False).std().to_numpy().flatten()
mag_sem = df[['mag','stim']].groupby(['stim'], sort=False).sem().to_numpy().flatten()
df[['mag','stim']].groupby(['stim'], sort=False).describe()


# In[18]:


plt.figure(figsize=(9,9))

plt.errorbar(np.arange(nori), mag_median, yerr=mag_sem, zorder=100, 
             color='silver', linewidth=3, capsize=3, capthick=3);
sns.stripplot(x="stim", y="mag", data=df);

# plt.ylim(-5,5);
plt.ylabel('adaptation magnitude');
plt.xlabel('stimulus identity');


# ### add dfof_ad thresholding

# In[19]:


sns.histplot(data=df, x='ad');


# In[20]:


sum(df.ad < 0) / len(df) # percent of negative dfof_ad


# In[21]:


plt.figure(figsize=(9,6))
ax = sns.scatterplot(data=df, x="ad", y="mag", hue="stim")
plt.xlim([-0.005,0.001]);
ax.set(title = 'adp mag changes with dfof_ad');


# In[22]:


# df2 = df.sort_values(by=['ad'])

# def f(win):
#     plt.figure(figsize=(15,5))
#     plt.plot(df2['ad'], df2['mag'].rolling(win, min_periods=1).mean(), alpha=0.7)
#     plt.xlim([-0.002,0.002])
#     plt.xlabel('dfof_ad')
#     plt.ylabel('adaptation mag rolling mean')
#     plt.title('adp mag rolling mean change with dfof_ad of cells')
#     plt.show()

# interactive_plot = interactive(f, win=(2, 20))
# output = interactive_plot.children[-1]
# output.layout.height = '350px'
# interactive_plot


# In[23]:


# df2 = df.sort_values(by=['ad'])

# def f(win):
#     plt.figure(figsize=(15,5))
#     plt.plot(df2['ad'], df2['mag'].rolling(win, min_periods=1).std(), alpha=0.7)
#     plt.xlim([-0.002,0.002])
#     plt.xlabel('dfof_ad')
#     plt.ylabel('adaptation mag rolling std')
#     plt.title('adp mag rolling std change with dfof_ad of cells')
#     plt.show()

# interactive_plot = interactive(f, win=(2, 20))
# output = interactive_plot.children[-1]
# output.layout.height = '350px'
# interactive_plot


# cell selection: vis_ad only, with dfof_ad thresholding

# In[24]:


dfof_threshold = 0.002
df_th = df.copy()
df_th.loc[df_th[np.abs(df.ad) < dfof_threshold].index.to_numpy(),'mag'] = np.nan
sns.histplot(data=df_th, x='mag');


# In[25]:


df_th[['mag','ad','tg','stim']].groupby(['stim'], sort=False).describe()


# In[26]:


mag_mean = df_th[['mag','stim']].groupby(['stim'], sort=False).mean().to_numpy().flatten()
mag_median = df_th[['mag','stim']].groupby(['stim'], sort=False).median().to_numpy().flatten()
mag_std = df_th[['mag','stim']].groupby(['stim'], sort=False).std().to_numpy().flatten()
mag_sem = df_th[['mag','stim']].groupby(['stim'], sort=False).sem().to_numpy().flatten()

plt.figure(figsize=(9,9))
plt.errorbar(np.arange(nori), mag_median, yerr=mag_sem, zorder=100, 
             color='silver', linewidth=3, capsize=3, capthick=3);
sns.stripplot(x="stim", y="mag", data=df_th);
plt.ylabel('adaptation magnitude');
plt.xlabel('stimulus identity');


# In[27]:


plt.figure(figsize=(9,6))
plt.errorbar(np.arange(nori)-0.03, mag_median, yerr=mag_sem, zorder=100, 
             color='gray', linewidth=3, capsize=3, capthick=3);
plt.errorbar(np.arange(nori)+0.03, mag_mean, yerr=mag_sem, zorder=100, 
             color='blue', linewidth=3, capsize=3, capthick=3);
plt.gca().set_ylim(top=0);
plt.legend(['median', 'mean'])
plt.ylabel('adaptation magnitude');
plt.xlabel('stimulus identity');


# In[40]:


mag_sort_id = mag_median.argsort()
mag_median_sorted = mag_median[mag_sort_id[::-1]]
mag_sem_sorted_med = mag_sem[mag_sort_id[::-1]]

mag_sort_id = mag_mean.argsort()
mag_mean_sorted = mag_mean[mag_sort_id[::-1]]
mag_sem_sorted_mean = mag_sem[mag_sort_id[::-1]]

mag_median_sorted[0], mag_median_sorted[-1], mag_mean_sorted[0], mag_mean_sorted[-1]


# In[42]:


plt.figure(figsize=(15,9))
plt.errorbar(np.arange(nori), mag_median_sorted, yerr=mag_sem_sorted_med, zorder=2, 
             color='silver', linewidth=3, capsize=3, capthick=3, alpha=0.5);
plt.errorbar(np.arange(nori), mag_mean_sorted, yerr=mag_sem_sorted_mean, zorder=1, 
             color='blue', linewidth=3, capsize=3, capthick=3, alpha=0.8);
plt.legend(['median', 'mean'])
plt.ylabel('adaptation magnitude');
plt.xlabel('stimulus identity sorted by adp mag');
# plt.ylim(-2.3,0.1);
plt.title('real data');


# ### compare w shuffle

# In[39]:


dfof_threshold = 0.0005
adp_threshold = 20 # 2

adp_shuffle = adp_mag.copy()
# print(adp_mag.shape) # ncell x nstim. should shuffle along rows (keep cell order, permute stim order)
adp_shuffle[np.abs(adp_shuffle) > adp_threshold] = np.nan

[np.random.shuffle(x) for x in adp_shuffle]; # shuffle every row

mag = adp_shuffle.flatten('F')

ad = dfof_ad.flatten('F')
mag[np.abs(ad) < dfof_threshold] = np.nan

stim = [np.arange(nori)] * adp_mag.shape[0]
stim_flat = np.sort([item for sublist in stim for item in sublist])

df_shuffle = pd.DataFrame({'stim':stim_flat, 'mag':mag, 'abs_mag':np.abs(mag)})
# df_shuffle.tail(5)

mag_shuffle_mean = df_shuffle[['mag','stim']].groupby(['stim'], sort=False).mean().to_numpy().flatten()
mag_shuffle_median = df_shuffle[['mag','stim']].groupby(['stim'], sort=False).median().to_numpy().flatten()
mag_shuffle_std = df_shuffle[['mag','stim']].groupby(['stim'], sort=False).std().to_numpy().flatten()
mag_shuffle_sem = df_shuffle[['mag','stim']].groupby(['stim'], sort=False).sem().to_numpy().flatten()
# df_shuffle[['mag','stim']].groupby(['stim'], sort=False).describe()

mag_shuffle_sort_id = mag_shuffle_median.argsort()
mag_shuffle_median_sorted = mag_shuffle_median[mag_shuffle_sort_id[::-1]]
mag_shuffle_sem_sorted_med = mag_shuffle_sem[mag_shuffle_sort_id[::-1]]

mag_shuffle_sort_id = mag_shuffle_mean.argsort()
mag_shuffle_mean_sorted = mag_shuffle_mean[mag_shuffle_sort_id[::-1]]
mag_shuffle_sem_sorted_mean = mag_shuffle_sem[mag_shuffle_sort_id[::-1]]
# mag_shuffle_median_sorted[0], mag_shuffle_median_sorted[-1], mag_shuffle_mean_sorted[0], mag_shuffle_mean_sorted[-1]

# plt.figure(figsize=(15,9))
# plt.errorbar(np.arange(nori), mag_shuffle_median_sorted, yerr=mag_shuffle_sem_sorted_med, zorder=2, 
#              color='silver', linewidth=3, capsize=3, capthick=3, alpha=0.5);
# plt.errorbar(np.arange(nori), mag_shuffle_mean_sorted, yerr=mag_shuffle_sem_sorted_mean, zorder=1, 
#              color='blue', linewidth=3, capsize=3, capthick=3, alpha=0.8);
# plt.legend(['median', 'mean'])
# plt.ylabel('adaptation mag_shufflenitude');
# plt.xlabel('fake stimulus identity sorted by adp mag');
# plt.title('shuffled stim');


# In[48]:


plt.figure(figsize=(15,9))
plt.errorbar(np.arange(nori), mag_shuffle_median_sorted, yerr=mag_shuffle_sem_sorted_med, zorder=2, 
             color='red', linewidth=3, capsize=3, capthick=3, alpha=0.2);
plt.errorbar(np.arange(nori), mag_median_sorted, yerr=mag_sem_sorted_med, zorder=1, 
             color='blue', linewidth=3, capsize=3, capthick=3, alpha=0.2);
plt.legend(['shuffle', 'real'])
plt.ylabel('adaptation mag median');
plt.xlabel('stimulus identity sorted by adp mag');
plt.gca().set_ylim(top=0.1)


# In[46]:


plt.figure(figsize=(15,9))
plt.errorbar(np.arange(nori), mag_shuffle_mean_sorted, yerr=mag_shuffle_sem_sorted_mean, zorder=2, 
             color='red', linewidth=3, capsize=3, capthick=3, alpha=0.2);
plt.errorbar(np.arange(nori), mag_mean_sorted, yerr=mag_sem_sorted_mean, zorder=1, 
             color='blue', linewidth=3, capsize=3, capthick=3, alpha=0.2);
plt.legend(['shuffle', 'real'])
plt.ylabel('adaptation mag mean');
plt.xlabel('stimulus identity sorted by adp mag');
plt.gca().set_ylim(top=0.2)


# ### dfof_ad & dfof_tg across stims

# In[288]:


df_th.loc[df_th[np.abs(df.ad) < dfof_threshold].index.to_numpy(),['ad','tg','mag']] = np.nan # threshold by dfof_ad
df_th.tail()


# In[307]:


stim_flat_cp = np.concatenate([df_th.stim, df_th.stim])
dfof_flat = np.concatenate([df_th.ad, df_th.tg])
adtg_flat = [['ad','tg']] * df_th.shape[0]
adtg_flat = np.sort([item for sublist in adtg_flat for item in sublist])
df_cp = pd.DataFrame({'stim':stim_flat_cp, 'dfof': dfof_flat, 'ad_or_tg': adtg_flat})
df_cp


# In[314]:


fig, ax = plt.subplots(figsize=(9,6))
ax = sns.barplot(x="stim", y="dfof", hue="ad_or_tg", data=df_cp, ci=68, palette="pastel") # mean +- sem


# ## Trace

# In[235]:


trace_mean = np.mean(trace, axis=0)
trace_std = np.std(trace, axis=0)
trace_sem = trace_std / np.sqrt(trace.shape[0])


# In[317]:


# fig, ax = plt.subplots(figsize=(20,15))
fig, ax = plt.subplots(figsize=(9,6))

for iori in np.arange(nori):
    ax.plot(np.arange(trace.shape[2]), trace_mean[iori,:].flatten(), color=sns.color_palette("bright")[iori])
    ax.fill_between(np.arange(trace.shape[2]), 
                    trace_mean[iori,:].flatten() + trace_sem[iori,:].flatten(), 
                    trace_mean[iori,:].flatten() - trace_sem[iori,:].flatten(),
                    color=sns.color_palette("bright")[iori], alpha=0.1)
plt.grid('minor')
plt.ylim(-0.02,0.10)
plt.xlabel('frame number')
plt.ylabel('dfof')
plt.show()


# ## Dim Reduct
# discard neurons with abs_adp > 2  
# then use full data, unadapted trials, or adapted trials to do PCA / TSNE

# In[39]:


df.describe()


# In[42]:


df = df[df.abs_mag <= 2]
df.describe()


# In[43]:


(764-726)/764 # percentage of discarded neurons


# In[44]:


plt.hist(df.abs_mag, bins = 100)
print("Median abs adaptation magnitude: ", df.abs_mag.median())

adapted_cells = df[df.abs_mag > df.abs_mag.median()]
less_adapted_cells = df[df.abs_mag <= df.abs_mag.median()]

plt.axvline(df.abs_mag.median(), color = 'red')


# In[45]:


df


# ### TSNE full data
# construct feature matrix

# In[89]:


temp = scipy.io.loadmat(os.path.join(dir_name, dir_sub, 'dim_red' + '.mat'))
temp.keys()


# In[64]:


df_adp = pd.DataFrame(temp['adp_vis'])
df_adp.columns = ['vis_ad', 'adp_abs']
df_adp.tail()


# In[55]:


plt.hist(df_adp.adp_abs, 100);


# In[66]:


df_adp.loc[df_adp.adp_abs > 5, 'adp_abs'] = np.nan
df_adp.loc[df_adp.vis_ad == 0, 'adp_abs'] = np.nan
df_adp


# In[73]:


cell_id = np.arange(len(df_adp))
df_adp['cell_id'] = cell_id
df_adp.drop(columns='vis_ad', inplace=True)
df_adp


# In[76]:


df = df_adp[~np.isnan(df_adp.adp_abs)].copy()
df


# In[77]:


plt.hist(df.adp_abs, bins = 100)
adp_median = np.median(df.adp_abs)
print("Median adaptation magnitude: ", adp_median)
plt.axvline(adp_median, color = 'red')


# In[95]:


adapt_cells = df[df.adp_abs > adp_median].cell_id.values
robust_cells = df[df.adp_abs <= adp_median].cell_id.values
robust_cells


# In[109]:


feature_ad = pd.DataFrame(temp['feature_ad']) # ntrial x ncell of adapter response
feature_tg = pd.DataFrame(temp['feature_tg']) # ntrial x ncell of target response
feature_full = pd.concat([feature_ad, feature_tg], axis=0) # full response

dim_red_ad = pd.DataFrame(temp['dim_red_ad']) # ntrial x ncell of adapter response
dim_red_tg = pd.DataFrame(temp['dim_red_tg']) # ntrial x ncell of target response
dim_red_full = pd.concat([dim_red_ad, dim_red_tg], axis=0) # full response


# In[110]:


adapt_cells_features = feature_full.loc[:, adapt_cells]
robust_cells_features = feature_full.loc[:, robust_cells]

adapt_cells_df = dim_red_full.loc[:, adapt_cells]
robust_cells_df = dim_red_full.loc[:, robust_cells]


# In[117]:


dim_red_full.rename(columns={0:'img_id'}, inplace=True)
dim_red_full.rename(columns={1:'rep_num'}, inplace=True)
dim_red_full[['img_id', 'rep_num']]


# In[118]:


adapt_tsne = TSNE(n_components=2, init='random',
                  random_state=0, perplexity=30, n_iter = 3000)


adapt_tsne_fit = adapt_tsne.fit_transform(adapt_cells_features)

adapt_tsne_df = pd.DataFrame(data = adapt_tsne_fit, columns = ['comp1', 'comp2']) #put the labels (image index and repeat number) back for plotting

adapt_tsne_df = pd.concat([adapt_tsne_df, dim_red_full[['img_id', 'rep_num']].reset_index()], axis = 1)


# In[120]:


plt.figure(figsize = (8,8))
plt.title('2 component t-SNE', fontsize = 20)
plt.scatter(adapt_tsne_df.comp1, adapt_tsne_df.comp2, c = adapt_tsne_df.img_id, cmap = 'Dark2', alpha = 0.7)

print(adapt_tsne.kl_divergence_)


# In[121]:


robust_tsne = TSNE(n_components=2, init='random',
                   random_state=0, perplexity=30, n_iter = 3000)
robust_tsne_fit = robust_tsne.fit_transform(robust_cells_features)

robust_tsne_df = pd.DataFrame(data =robust_tsne_fit, columns = ['comp1', 'comp2'])
robust_tsne_df  = pd.concat([robust_tsne_df, dim_red_full[['img_id', 'rep_num']].reset_index()], axis = 1)

plt.figure(figsize = (8,8))
plt.title('2 component t-SNE', fontsize = 20)
plt.scatter(robust_tsne_df.comp1, robust_tsne_df.comp2, c = robust_tsne_df.img_id, cmap = 'Dark2', alpha = 0.7)

print(robust_tsne.kl_divergence_)


# ### PCA, full data

# In[123]:


above_pca = PCA(n_components = 3)
above_principal_pcs = above_pca.fit_transform(adapt_cells_features)

above_principal_df = pd.DataFrame(data = above_principal_pcs, columns = ['pc1', 'pc2', 'pc3'])
above_principal_df = pd.concat([above_principal_df, dim_red_full[['img_id', 'rep_num']].reset_index()], axis = 1)

fig = plt.figure(figsize = (8,8))
ax = plt.axes(projection='3d')
ax.scatter3D(above_principal_df.pc1, above_principal_df.pc2, above_principal_df.pc3, c = above_principal_df.img_id, cmap = 'Dark2' )
#ax.view_init(60, 30) 


# In[122]:


below_pca = PCA(n_components = 3)
below_principal_pcs = below_pca.fit_transform(robust_cells_features)

below_principal_df = pd.DataFrame(data = below_principal_pcs, columns = ['pc1', 'pc2', 'pc3'])
below_principal_df = pd.concat([below_principal_df, dim_red_full[['img_id', 'rep_num']].reset_index()], axis = 1)

fig = plt.figure(figsize = (8,8))
ax = plt.axes(projection='3d')
ax.scatter3D(below_principal_df.pc1, below_principal_df.pc2, below_principal_df.pc3, c = below_principal_df.img_id, cmap = 'Dark2' )
#ax.view_init(60, 30) 


# In[127]:


above_pca = PCA()
below_pca = PCA()

above_pcs = above_pca.fit_transform(adapt_cells_features)
below_pcs = below_pca.fit_transform(robust_cells_features)

plt.plot(np.cumsum(above_pca.explained_variance_ratio_))
plt.plot(np.cumsum(below_pca.explained_variance_ratio_))
plt.legend(['adapt', 'robust']);


# ## Misc
# ### Adaptation increases variability for natural stim too
# for target after adaptation, fano factor is higher than adapter  
# cell selection: vis_ad only, with thresholding of dfof_ad

# In[270]:


adp_fano_tg = dfof_tg_std / dfof_tg
adp_fano_ad = dfof_ad_std / dfof_ad
adp_fano = (adp_fano_tg - adp_fano_ad) / (adp_fano_tg + adp_fano_ad)


# In[274]:


adp_fano_tg = dfof_tg_std / dfof_tg
adp_fano_ad = dfof_ad_std / dfof_ad
adp_fano = (adp_fano_tg - adp_fano_ad) / (adp_fano_tg + adp_fano_ad) # range [-1,1], meaning fano factor de/increase after adp

df_fano = df_th.copy()
df_fano['fano'] = adp_fano.flatten('F')
df_fano.tail()


# In[273]:


df_fano[['fano','stim']].groupby(['stim'], sort=False).describe()


# In[278]:


fano_mean = df_fano[['fano','stim']].groupby(['stim'], sort=False).mean().to_numpy().flatten()
fano_median = df_fano[['fano','stim']].groupby(['stim'], sort=False).median().to_numpy().flatten()
fano_std = df_fano[['fano','stim']].groupby(['stim'], sort=False).std().to_numpy().flatten()
fano_sem = df_fano[['fano','stim']].groupby(['stim'], sort=False).sem().to_numpy().flatten()

plt.figure(figsize=(9,9))
plt.errorbar(np.arange(nori), fano_median, yerr=fano_sem, zorder=100, 
             color='silver', linewidth=3, capsize=3, capthick=3);
sns.stripplot(x="stim", y="fano", data=df_fano);
plt.ylabel('delta fano factor after adaptation');
plt.xlabel('stimulus identity');
plt.ylim(-2,2);


# ### stim presentation order

# In[12]:


dir_behav = 'Z:/All_Staff/home/lan/Data/2P_images/i1338/210616'
behav = scipy.io.loadmat(os.path.join(dir_behav, 'input_behav.mat'))


# In[24]:


flat_list = [item for sublist in behav['stim1'] for item in sublist]
flat_list = [item for sublist in flat_list for item in sublist]
flat_list = [item for sublist in flat_list for item in sublist]
stim1 = np.array(flat_list)

# flat_list = [item for sublist in behav['stim2'] for item in sublist]
# flat_list = [item for sublist in flat_list for item in sublist]
# flat_list = [item for sublist in flat_list for item in sublist]
# stim2 = np.array(flat_list)
# stim2 == stim1


# In[31]:


values, counts = np.unique(stim1, return_counts=True)
counts


# In[ ]:




