#!/usr/bin/env python
# coding: utf-8

# In[1]:


import os
import numpy as np
import pandas as pd
import scipy.io

import seaborn as sns
import matplotlib.pyplot as plt
# from ipywidgets import interactive
get_ipython().run_line_magic('matplotlib', 'inline')


# ## Load data

# In[2]:


d = {'mouse': [1329, 1329], 'date': ['201209', '201209_caiman'], 'area': ['V1', 'V1']}
meta = pd.DataFrame(data=d)
meta


# In[3]:


nset = len(meta.index); ncell = []; nori = 8; nisi = 3; nframe_trial = 207
# dir_name = 'C:\\Users\\lan\\Documents\\repos\\inter\\mat\\'
root_folder = 'C://Users/ll357/Documents/inter/'
dir_name = root_folder + 'mat/'

vis_ad = np.empty([0,1]); vis_tg = np.empty([0,1]); well_fit = np.empty([0,1])
ori_pref = np.empty([0,nisi]); fit_param = np.empty([0,7,nisi])
dfof_ad = np.empty([0,1]); dfof_tg = np.empty([0,nori,nisi])
dfof_ad_std = np.empty([0,1]); dfof_tg_std = np.empty([0,nori,nisi])
trace = np.empty([0,nori,nisi,nframe_trial])

for iset in np.arange(nset):
    dir_sub = str(meta.area[iset]) + '_i' + str(meta.mouse[iset]) + '_' + str(meta.date[iset])

    cell_prop = scipy.io.loadmat(os.path.join(dir_name, dir_sub, 'cell_property_loose' + '.mat'))
    dfof = scipy.io.loadmat(os.path.join(dir_name, dir_sub, 'dfof' + '.mat'))
    trace_align = scipy.io.loadmat(os.path.join(dir_name, dir_sub, 'trace_aligned' + '.mat'))
    fit_tuning = scipy.io.loadmat(os.path.join(dir_name, dir_sub, 'fit_tuning_isi3' + '.mat'))
    
    ncell.append(len(cell_prop['vis_cell_ad']))    
    vis_ad = np.concatenate((vis_ad, cell_prop['vis_cell_ad']), axis=0)
    vis_tg = np.concatenate((vis_tg, cell_prop['vis_cell_noad_tg']), axis=0)
#     well_fit = np.concatenate((well_fit, cell_prop['well_fit_cell']), axis=0)
    
    ori_pref = np.concatenate((ori_pref, cell_prop['ori_pref']), axis=0)
    fit_param = np.concatenate((fit_param, fit_tuning['fit_param']), axis=0)

    dfof_ad = np.concatenate((dfof_ad, dfof['dfof_ad']), axis=0)
    dfof_ad_std = np.concatenate((dfof_ad_std, dfof['dfof_ad_std']), axis=0)
    dfof_tg = np.concatenate((dfof_tg, dfof['dfof_tg']), axis=0)
    dfof_tg_std = np.concatenate((dfof_tg_std, dfof['dfof_tg_std']), axis=0)
    
    trace_flat = np.empty([ncell[iset],nori,nisi,nframe_trial]);
    for icell in np.arange(ncell[iset]):
        for iori in np.arange(nori):
            for iisi in np.arange(nisi):
                trace_flat[icell][iori][iisi][:] = trace_align['trace_avg'][icell][iori][iisi].flatten()
    trace = np.vstack((trace,trace_flat))

ncell, vis_ad.shape, vis_tg.shape#, well_fit.shape, ori_pref.shape, fit_param.shape, dfof_ad.shape, dfof_tg.shape, trace.shape


# In[4]:


meta['ncell'] = ncell
meta.loc[1, 'area'] = 'V1_caiman'
meta


# In[5]:


mouse_cell = [item for item, count in zip(meta.mouse, meta.ncell) for i in range(count)]
area_cell = [item for item, count in zip(meta.area, meta.ncell) for i in range(count)]
meta_cell = pd.DataFrame({'mouse': mouse_cell, 'area': area_cell})
meta_cell


# In[6]:


meta_cell.area.value_counts()


# ## Adaptation magnitude

# adaptation mag = response to target ori==0 with adapter / response to adapter - 1  
# cell selection: vis_ad only, no dfof_ad thresholding

# In[7]:


adp_mag = dfof_tg[:,0,[1,2]] / dfof_ad - 1 # isi=750 -> 250

meta_cell_750 = meta_cell.copy(); meta_cell_750['isi'] = 750
meta_cell_250 = meta_cell.copy(); meta_cell_250['isi'] = 250
meta_cell_isi = pd.concat([meta_cell_750, meta_cell_250], ignore_index=True)

df_adp_mag = meta_cell_isi.copy()
df_adp_mag['adp_mag'] = adp_mag.flatten('F') # flatten by concat columns
df_adp_mag['dfof_ad'] = np.concatenate((dfof_ad, dfof_ad), axis=0)

df_adp_mag['vis_ad'] = np.concatenate((vis_ad, vis_ad), axis=0)
df_adp_mag = df_adp_mag[ df_adp_mag['vis_ad'] == 1 ]
df_adp_mag#.reset_index()


# In[8]:


df_adp_mag[df_adp_mag.area == 'V1_caiman'].vis_ad.describe() # 151 caiman cell are vis driven by ad


# In[9]:


df_adp_mag[['adp_mag','area','isi']].groupby(['area','isi'], sort=False).describe()


# In[10]:


sns.set_theme(style="whitegrid")
sns.set_context("talk", font_scale = 0.9)

plt.figure(figsize=(8,8))
ax = sns.violinplot(data=df_adp_mag, x="area", y="adp_mag", hue="isi",
               split=True, inner="quart", palette="Set3")
sns.despine(left=True)
ax.set(title = 'adaptation magnitude before dfof_ad thresholding')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# ### add dfof_ad thresholding

# In[11]:


g = sns.FacetGrid(df_adp_mag, col="area",  row="isi")
g.map(sns.scatterplot, "dfof_ad", "adp_mag")
plt.xlim([0,0.05]);


# In[12]:


# plt.figure(figsize=(9,6))
# ax = sns.scatterplot(data=df_adp_mag[df_adp_mag.area == 'V1_caiman'], x="dfof_ad", y="adp_mag", hue="area", style="isi")
# # plt.xlim([0,0.2]);
# plt.xlim([0,0.05]);
# ax.set(title = 'adp mag changes with dfof_ad');


# In[13]:


# df = df_adp_mag.sort_values(by=['dfof_ad'])
# df1 = df[df.isi == 750]
# df2 = df[df.isi == 250]

# def f(win):
#     plt.figure(figsize=(15,5))
#     plt.plot(df1.dfof_ad, df1['adp_mag'].rolling(win, min_periods=1).mean(), alpha=0.7)
#     plt.plot(df2.dfof_ad, df2['adp_mag'].rolling(win, min_periods=1).mean(), alpha=0.7)
#     plt.legend(['isi = 750', 'isi = 250'])
#     plt.xlim([0,0.1])
#     plt.xlabel('dfof_ad')
#     plt.ylabel('adaptation mag rolling mean')
#     plt.title('adp mag rolling mean change with dfof_ad of cells')
#     plt.show()

# interactive_plot = interactive(f, win=(2, 20))
# output = interactive_plot.children[-1]
# output.layout.height = '350px'
# interactive_plot


# In[14]:


# def f(win):
#     plt.figure(figsize=(15,5))
#     plt.plot(df1.dfof_ad, df1['adp_mag'].rolling(win, min_periods=1).std(), alpha=0.7)
#     plt.plot(df2.dfof_ad, df2['adp_mag'].rolling(win, min_periods=1).std(), alpha=0.7)
#     plt.legend(['isi = 750', 'isi = 250'])
#     plt.xlim([0,0.1])
#     plt.xlabel('dfof_ad')
#     plt.ylabel('adaptation mag rolling std')
#     plt.title('adp mag rolling std change with dfof_ad of cells')
#     plt.show()

# interactive_plot = interactive(f, win=(2, 20))
# output = interactive_plot.children[-1]
# output.layout.height = '350px'
# interactive_plot


# In[15]:


plt.figure(figsize=(15,10));
plt.subplot(221)
sns.histplot(data=df_adp_mag[df_adp_mag.isi == 750], x="dfof_ad", hue="area", kde=True)
plt.subplot(222)
sns.histplot(data=df_adp_mag[df_adp_mag.isi == 750], x="adp_mag", hue="area", kde=True)
plt.subplot(223)
sns.histplot(data=df_adp_mag[df_adp_mag.isi == 250], x="dfof_ad", hue="area", kde=True)
plt.subplot(224)
sns.histplot(data=df_adp_mag[df_adp_mag.isi == 250], x="adp_mag", hue="area", kde=True)


# cell selection: vis_ad only, with dfof_ad thresholding

# In[105]:


dfof_threshold_manual = 0.03 # 0.03, 0.05
dfof_threshold_caiman = 0.02 # 0.02

df_adp_mag_thres = df_adp_mag[((np.abs(df_adp_mag.dfof_ad) >= dfof_threshold_manual) & (df_adp_mag.area == 'V1'))
                             | ((np.abs(df_adp_mag.dfof_ad) >= dfof_threshold_caiman) & (df_adp_mag.area == 'V1_caiman'))]
# df_adp_mag_thres.reset_index();
df_adp_mag_thres


# In[106]:


df_adp_mag_thres[['adp_mag','area','isi']].groupby(['area','isi'], sort=False).describe()


# In[107]:


plt.figure(figsize=(8,8))
ax = sns.violinplot(data=df_adp_mag_thres, x="area", y="adp_mag", hue="isi",
               split=True, inner="quart", palette="Set3")
sns.despine(left=True)
ax.set(ylabel = 'adaptation index', title = 'adaptation magnitude after dfof_ad thresholding');
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# In[108]:


plt.figure(figsize=(8,8))

ax = sns.boxplot(data=df_adp_mag_thres, x="area", y="adp_mag", hue="isi") #, boxprops=dict(alpha=.3))
for patch in ax.artists:
    r, g, b, a = patch.get_facecolor()
    patch.set_facecolor((r, g, b, .3))
ax = sns.swarmplot(data=df_adp_mag_thres, x="area", y="adp_mag", hue="isi", dodge=True)

sns.despine(left=True)
ax.set(ylabel = 'adaptation index', title = 'adaptation magnitude after dfof_ad thresholding');
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# In[109]:


sns.scatterplot(data=df_adp_mag_thres[df_adp_mag_thres.isi==250], x="dfof_ad", y="adp_mag", hue="area")
plt.xlim(0,0.1)


# In[110]:


# t = df_adp_mag_thres.copy()
# tt = t[['adp_mag','area','isi']].groupby(['area','isi'], sort=False).median().reset_index()

ax = sns.pointplot(data=df_adp_mag_thres, x="area", y="adp_mag", hue="isi", ci=68)
sns.despine()
ax.set(ylabel = 'adaptation index')
plt.ylim([-1,0])
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# ## Trace

# ### trace using all cells

# In[111]:


trace_mean = np.mean(trace, axis=0)
trace_std = np.std(trace, axis=0)
trace_sem = trace_std / np.sqrt(trace.shape[0])

# compare btw isi
trace_mean = np.mean(trace_mean, axis=0)
trace_std = np.std(trace_std, axis=0)
trace_sem = trace_std / np.sqrt(trace.shape[0])
trace_mean.shape


# In[80]:


fig, ax = plt.subplots(figsize=(9,6))

for iisi in np.arange(nisi):
    ax.plot(np.arange(trace.shape[-1]), trace_mean[iisi,:].flatten(), 
            color=sns.color_palette("bright")[iisi])
    
    ax.fill_between(np.arange(trace.shape[-1]), 
                    trace_mean[iisi,:].flatten() + trace_sem[iisi,:].flatten(), 
                    trace_mean[iisi,:].flatten() - trace_sem[iisi,:].flatten(),
                    color=sns.color_palette("bright")[iisi], alpha=0.1)
plt.grid('minor')
plt.xlim(0,70)
plt.xlabel('frame number')
plt.ylabel('dfof')
plt.legend(['inf', '750', '250'])
plt.show()


# In[97]:


trace_mean = np.mean(trace, axis=0)
trace_std = np.std(trace, axis=0)
trace_sem = trace_std / np.sqrt(trace.shape[0])

# visliz only isi=250
trace_mean = trace_mean[:,-1,:]
# trace_std = trace_std[:,-1,:]
trace_sem = trace_sem[:,-1,:]
trace_mean.shape


# In[82]:


fig, ax = plt.subplots(figsize=(9,6))

for iori in np.arange(nori):
    ax.plot(np.arange(trace.shape[-1]), trace_mean[iori,:].flatten(), color=sns.color_palette("bright")[iori])
    ax.fill_between(np.arange(trace.shape[-1]), 
                    trace_mean[iori,:].flatten() + trace_sem[iori,:].flatten(), 
                    trace_mean[iori,:].flatten() - trace_sem[iori,:].flatten(),
                    color=sns.color_palette("bright")[iori], alpha=0.1)
plt.grid('minor')
plt.xlim(0,70)
# plt.ylim(-0.02,0.10)
plt.xlabel('frame number')
plt.ylabel('dfof')
plt.show()


# ### trace using vis driven & post thres dfof_ad

# In[112]:


trace.shape, ncell, vis_ad.shape


# In[115]:


vis_id = df_adp_mag_thres[df_adp_mag_thres.isi==750].index.values
vis_id_manual = vis_id[vis_id < ncell[0]]
vis_id_caiman = vis_id[vis_id >= ncell[0]]
print(vis_id_manual)
print(vis_id_caiman)
len(vis_id_manual), len(vis_id_caiman)


# In[116]:


print(ncell[0])
vis_id_caiman_0 = vis_id_caiman - ncell[0]
vis_id_caiman_0


# In[102]:


# import pickle

# f = open('vis_cell_manual_vs_caiman.pckl', 'wb')
# pickle.dump([vis_id_manual, vis_id_caiman_0], f)
# f.close()

# # f = open('vis_cell_caiman_vs_manual.pckl', 'rb')
# # vis_id_manual, vis_id_caiman_0 = pickle.load(f)
# # f.close()


# In[117]:


# select trace from vis-driven caiman / manual cells

vis_id = df_adp_mag_thres[df_adp_mag_thres.isi==750].index.values
vis_id_manual = vis_id[vis_id < ncell[0]]
vis_id_caiman = vis_id[vis_id >= ncell[0]]
# print(vis_id_manual)
# print(vis_id_caiman)

trace_manual = trace[vis_id_manual, :,:,:]
trace_caiman = trace[vis_id_caiman, :,:,:]
trace_manual.shape, trace_caiman.shape


# In[120]:


fig, ax = plt.subplots(1,2,figsize=(12,6)) # , sharey=True
ax_id = 0

for trace_now in [trace_manual, trace_caiman]:
    trace_mean = np.mean(np.mean(trace_now, axis=0), axis=0)
    trace_flatter = trace_now.reshape([-1, trace_now.shape[-2], trace_now.shape[-1]]) # flatten ncell x nori
    trace_std = np.std(trace_flatter, axis=0)
    trace_sem = trace_std / np.sqrt(trace_flatter.shape[0])

    for iisi in np.arange(nisi):
        ax[ax_id].plot(np.arange(trace_now.shape[-1]), trace_mean[iisi,:].flatten(), 
                color=sns.color_palette("bright")[iisi])

        ax[ax_id].fill_between(np.arange(trace_now.shape[-1]), 
                        trace_mean[iisi,:].flatten() + trace_sem[iisi,:].flatten(), 
                        trace_mean[iisi,:].flatten() - trace_sem[iisi,:].flatten(),
                        color=sns.color_palette("bright")[iisi], alpha=0.1)

    ax[ax_id].set_xlim(0,70)
    ax[0].set_ylim(-0.01, 0.13)
    ax[1].set_ylim(-0.01, 0.13)
    ax[ax_id].set_xlabel('frame number')
    ax[ax_id].set_ylabel('dfof')
    ax[0].set_title('manual trace')
    ax[1].set_title('caiman trace')
    plt.legend(['inf', '750', '250'])
    
    ax_id = ax_id + 1


# ## Adaptation tuning bias

# In[53]:


tt = ori_pref.copy()
tt[tt > 90] = np.abs(tt[tt > 90] - 180)
tuning_bias = tt[:,[1,2]] - tt[:,[0]];

ori_pref_bin = tt[:,[0]];
ori_pref_bin[ori_pref_bin < 22.5] = 0; ori_pref_bin[ori_pref_bin > 67.5] = 90; 
ori_pref_bin[(ori_pref_bin >= 22.5) & (ori_pref_bin <= 67.5)] = 45; 


# ### run well fit func first

# In[55]:


df_adp_tune = meta_cell_isi.copy()
df_adp_tune['tuning_bias'] = tuning_bias.flatten('F')
df_adp_tune['ori_pref_bin'] = np.concatenate((ori_pref_bin, ori_pref_bin), axis=0)

df_adp_tune['vis_ad'] = np.concatenate((vis_ad, vis_ad), axis=0)
# df_adp_tune['well_fit'] = np.concatenate((well_fit, well_fit), axis=0)
df_adp_tune = df_adp_tune[ df_adp_tune['vis_ad'] == 1 ]
# df_adp_tune = df_adp_tune[ df_adp_tune['well_fit'] == 1 ]

b, c = df_adp_tune.iloc[0].copy(), df_adp_tune.iloc[1].copy() 
df_adp_tune.iloc[0], df_adp_tune.iloc[1] = c, b # swap row 0 & 1 to sort df.gb properly
df_adp_tune.reset_index()


# In[29]:


df_adp_tune[['tuning_bias','area','isi','ori_pref_bin']].groupby(['area','isi','ori_pref_bin'], sort=False).describe()


# In[30]:


sns.set_style("whitegrid")
plt.figure(figsize=(8,5))
ax = sns.violinplot(data=df_adp_tune[df_adp_tune.area == 'V1'], 
                    x="ori_pref_bin", y="tuning_bias", hue="isi", 
                    split=True, inner="quart", palette="Set3")
sns.despine(left=True)
ax.set(xlabel = '|pref - adapter|', ylabel = 'tuning bias', title = 'V1 tuning bias after adaptation')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# In[31]:


tt = df_adp_tune[['tuning_bias','area','isi','ori_pref_bin']].groupby(['area','isi','ori_pref_bin'], sort=False).median().reset_index()

ax = sns.lineplot(data=tt, x="ori_pref_bin", y="tuning_bias", hue="area", style="isi");
ax.set(xlabel = '|pref - adapter|', ylabel = 'tuning bias', title = 'median tuning bias after adaptation')
plt.xlim([0,90])
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# ## Misc
# ### Adaptation increases variability (?)
# for ori=0 target after adaptation, fano factor is higher than adapter  
# cell selection: vis_ad only, with thresholding of dfof_ad & dfof_tg0

# In[138]:


adp_fano_tg = dfof_tg_std[:,0,[1,2]] / dfof_tg[:,0,[1,2]]
adp_fano_ad = dfof_ad_std / dfof_ad
adp_fano = (adp_fano_tg - adp_fano_ad) / (adp_fano_tg + adp_fano_ad) # range [-1,1], meaning fano factor de/increase after adp
# adp_fano = (adp_fano_tg - adp_fano_ad) / (adp_fano_ad)

df_adp_fano = meta_cell_isi.copy()
df_adp_fano['adp_fano'] = adp_fano.flatten('F')
df_adp_fano['dfof_ad'] = np.concatenate((dfof_ad, dfof_ad), axis=0)
df_adp_fano['dfof_tg'] = np.concatenate((dfof_tg[:,0,1], dfof_tg[:,0,2]), axis=0)

df_adp_fano['vis_ad'] = np.concatenate((vis_ad, vis_ad), axis=0)
df_adp_fano = df_adp_fano[ df_adp_fano['vis_ad'] == 1 ]

df_adp_fano = df_adp_fano[(df_adp_fano.dfof_ad >= dfof_threshold) & (df_adp_fano.dfof_tg >= dfof_threshold)]
df_adp_fano.reset_index()


# In[139]:


df_adp_fano[['adp_fano','area','isi']].groupby(['area','isi'], sort=False).describe()


# In[140]:


sns.set_style("whitegrid")
plt.figure(figsize=(8,5))
ax = sns.violinplot(data=df_adp_fano, x="area", y="adp_fano", hue="isi",
               split=True, inner="quart", palette="Set3")
sns.despine(left=True)
ax.set(ylabel = 'fano factor change after adp', title = 'adaptation impacts variability')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# In[ ]:




