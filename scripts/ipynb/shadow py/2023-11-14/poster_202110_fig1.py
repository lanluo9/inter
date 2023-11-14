#!/usr/bin/env python
# coding: utf-8

# # Fig 1B

# In[33]:


import os
import numpy as np
import pandas as pd
import scipy.io

import seaborn as sns
sns.set_theme(style="whitegrid")
sns.set_context("talk", font_scale = 0.8)

import matplotlib.pyplot as plt
get_ipython().run_line_magic('matplotlib', 'inline')

# from ipywidgets import interactive


# In[2]:


from custom_functions import loadmat
help(loadmat)


# ## Load data

# In[3]:


pwd


# In[4]:


root_folder = 'C://Users/ll357/Documents/inter/'
# root_folder = 'C:/Users/lan/Documents/repos/inter/'
master_xls_file = 'mat/adp_dataset_master.xlsx'
meta = pd.read_excel(root_folder + master_xls_file, index_col=None)

# meta = meta[meta.seg == 'segmented'].reset_index()
meta = meta[meta.py_plot_ready == 'yes'].reset_index()

meta = meta[['mouse','date','area']]
meta.mouse = meta.mouse.astype(int)
meta.date = meta.date.astype(int)


# In[5]:


nset = len(meta.index); ncell = []; nori = 8; nisi = 3; nframe_trial = 223
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
    well_fit = np.concatenate((well_fit, cell_prop['well_fit_cell']), axis=0)
    
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

ncell, vis_ad.shape, vis_tg.shape, well_fit.shape, ori_pref.shape, fit_param.shape, dfof_ad.shape, dfof_tg.shape, trace.shape


# In[6]:


meta['ncell'] = ncell
# meta = meta.replace({'area' : { 'V1':1, 'LM':2, 'LI':3 }})
meta


# In[7]:


mouse_cell = [item for item, count in zip(meta.mouse, meta.ncell) for i in range(count)]
area_cell = [item for item, count in zip(meta.area, meta.ncell) for i in range(count)]
meta_cell = pd.DataFrame({'mouse': mouse_cell, 'area': area_cell})
meta_cell


# ## Adaptation magnitude

# adaptation mag = response to target ori==0 with adapter / response to adapter - 1  
# cell selection: vis_ad only, no dfof_ad thresholding

# In[8]:


adp_mag = dfof_tg[:,0,[1,2]] / dfof_ad - 1

meta_cell_750 = meta_cell.copy(); meta_cell_750['isi'] = 750
meta_cell_250 = meta_cell.copy(); meta_cell_250['isi'] = 250
meta_cell_isi = pd.concat([meta_cell_750, meta_cell_250], ignore_index=True)

df_adp_mag = meta_cell_isi.copy()
df_adp_mag['adp_mag'] = adp_mag.flatten('F')
df_adp_mag['dfof_ad'] = np.concatenate((dfof_ad, dfof_ad), axis=0)

df_adp_mag['vis_ad'] = np.concatenate((vis_ad, vis_ad), axis=0)
df_adp_mag = df_adp_mag[ df_adp_mag['vis_ad'] == 1 ]
df_adp_mag.reset_index()


# In[9]:


df_adp_mag[['adp_mag','area','isi']].groupby(['area','isi'], sort=False).describe()


# In[11]:


# sns.set_style("whitegrid")
# plt.figure(figsize=(8,5))
# ax = sns.violinplot(data=df_adp_mag, x="area", y="adp_mag", hue="isi",
#                split=True, inner="quart", palette="Set3")
# sns.despine(left=True)
# ax.set(title = 'adaptation magnitude before dfof_ad thresholding')
# plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# In[13]:


# sns.set_theme(style="whitegrid")
# sns.set_context("talk", font_scale = 0.8)

# plt.figure(figsize=(8,5))
# ax = sns.violinplot(data=df_adp_mag, x="area", y="adp_mag", hue="isi",
#                split=True, inner="quart", palette="Set3")
# sns.despine(left=True)
# ax.set(title = 'adaptation magnitude before dfof_ad thresholding')
# plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# ### add dfof_ad thresholding

# In[15]:


# plt.figure(figsize=(9,6))
# ax = sns.scatterplot(data=df_adp_mag, x="dfof_ad", y="adp_mag", hue="area", style="isi")
# # plt.xlim([0,0.8]);
# plt.xlim([0,0.1]);
# ax.set(title = 'adp mag changes with dfof_ad');


# In[17]:


df = df_adp_mag.sort_values(by=['dfof_ad'])
df1 = df[df.isi == 750]
df2 = df[df.isi == 250]

def f(win):
    plt.figure(figsize=(15,5))
    plt.plot(df1.dfof_ad, df1['adp_mag'].rolling(win, min_periods=1).mean(), alpha=0.7)
    plt.plot(df2.dfof_ad, df2['adp_mag'].rolling(win, min_periods=1).mean(), alpha=0.7)
    plt.legend(['isi = 750', 'isi = 250'])
    plt.xlim([0,0.1])
    plt.xlabel('dfof_ad')
    plt.ylabel('adaptation mag rolling mean')
    plt.title('adp mag rolling mean change with dfof_ad of cells')
    plt.show()

# interactive_plot = interactive(f, win=(2, 20))
# output = interactive_plot.children[-1]
# output.layout.height = '350px'
# interactive_plot


# In[18]:


def f(win):
    plt.figure(figsize=(15,5))
    plt.plot(df1.dfof_ad, df1['adp_mag'].rolling(win, min_periods=1).std(), alpha=0.7)
    plt.plot(df2.dfof_ad, df2['adp_mag'].rolling(win, min_periods=1).std(), alpha=0.7)
    plt.legend(['isi = 750', 'isi = 250'])
    plt.xlim([0,0.1])
    plt.xlabel('dfof_ad')
    plt.ylabel('adaptation mag rolling std')
    plt.title('adp mag rolling std change with dfof_ad of cells')
    plt.show()

# interactive_plot = interactive(f, win=(2, 20))
# output = interactive_plot.children[-1]
# output.layout.height = '350px'
# interactive_plot


# cell selection: vis_ad only, with dfof_ad thresholding

# In[19]:


dfof_threshold = 0.03
df_adp_mag_thres = df_adp_mag[df_adp_mag.dfof_ad >= dfof_threshold]
df_adp_mag_thres.reset_index()


# In[20]:


df_adp_mag_thres[['adp_mag','area','isi']].groupby(['area','isi'], sort=False).describe()


# In[21]:


# sns.set_style("whitegrid")
# plt.figure(figsize=(10,6))
# ax = sns.violinplot(data=df_adp_mag_thres, x="area", y="adp_mag", hue="isi",
#                split=True, inner="quart", palette="Set3")
# sns.despine(left=True)
# ax.set(ylabel = 'adaptation index', title = 'adaptation magnitude after dfof_ad thresholding');
# plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# In[48]:


adp_v1 = df_adp_mag_thres[(df_adp_mag_thres.isi==250) & (df_adp_mag_thres.area=='V1')].adp_mag.values
adp_lm = df_adp_mag_thres[(df_adp_mag_thres.isi==250) & (df_adp_mag_thres.area=='LM')].adp_mag.values
adp_li = df_adp_mag_thres[(df_adp_mag_thres.isi==250) & (df_adp_mag_thres.area=='LI')].adp_mag.values


# In[51]:


from scipy.stats import ttest_ind

t, p = ttest_ind(adp_v1, adp_lm, equal_var=False)
print("ttest_ind:            t = %g  p = %g" % (t, p))

t, p = ttest_ind(adp_li, adp_lm, equal_var=False)
print("ttest_ind:            t = %g  p = %g" % (t, p))


# ## fig 1B

# In[197]:


sns.set_style("ticks")
ax = sns.violinplot(x="area", y="adp_mag", col="isi",
                data=df_adp_mag_thres[df_adp_mag_thres.isi==250], kind="violin",
                height=4, aspect=1, palette='pastel');

# statistical annotation
x1, x2 = 0, 0.98
y, h, col = np.max(df_adp_mag_thres[df_adp_mag_thres.isi==250].adp_mag.values) + 0.4, 0.3, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, "****", ha='center', va='bottom', color=col, fontsize=16)

x1, x2 = 1.02, 2
y, h, col = np.max(df_adp_mag_thres[df_adp_mag_thres.isi==250].adp_mag.values) + 0.1, 0.3, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, "****", ha='center', va='bottom', color=col, fontsize=16)

x1, x2 = 0, 2
y, h, col = np.max(df_adp_mag_thres[df_adp_mag_thres.isi==250].adp_mag.values) + 1, 0.3, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, "****", ha='center', va='bottom', color=col, fontsize=16)

plt.xlabel('')
plt.ylabel('adaptation index', fontsize=18)
ax.set_xticklabels(['V1', 'LM', 'LI'], fontsize=18)
ax.set_yticks(np.arange(-3,3.1))
ax.set_ylim((-3,3))

ax.tick_params(axis="x", labelsize=16)
ax.tick_params(axis="y", labelsize=14)
sns.despine(ax=ax)

# figtitle = 'Adaptation magnitude increases along ventral stream'
# ax.fig.suptitle(figtitle, size=16)
# ax.fig.subplots_adjust(top=.8)
# ax.savefig(root_folder + 'plot/poster/' + figtitle + '.pdf', format='pdf')


# # Fig 1C

# In[120]:


import pandas as pd
import numpy as np
import scipy
from scipy.io import loadmat
from sklearn.linear_model import LinearRegression
import seaborn as sns
sns.set_style("whitegrid")
import matplotlib.pyplot as plt
from astropy.stats import circvar
from astropy import units as u
import os
# from ipywidgets import interactive


# In[121]:


def replace_area_name(df, area_key):
# df = dataframe containing an area array column
# area_key = string

    df[area_key] = df[area_key].replace(1, 'V1')
    df[area_key] = df[area_key].replace(2, 'LM')
    df[area_key] = df[area_key].replace(3, 'LI')
    
    return df


# In[122]:


mode = 'hubel'

if mode == 'nuke':
    root_path = "C:/Users/lan/Documents/repos/inter"
elif mode == 'hubel':
    root_path = 'C:/Users/ll357/Documents/inter/'
elif mode == 'linux':
    root_path = '/home/ll357@dhe.duke.edu/inter'

mode, root_path


# ## response amplitude by area

# In[123]:


file_path = os.path.join(root_path, "plot/CV SNR OSI R2 ori_perc by area - why HVA lack well fit", "amp_area.mat").replace("\\","/")
temp = loadmat(file_path)

df = pd.DataFrame(temp['resp_ad'], columns=['resp_ad'])
df['area_ad'] = temp['area_ad']
df['area_ad'] = df['area_ad'].replace(1, 'V1')
df['area_ad'] = df['area_ad'].replace(2, 'LM')
df['area_ad'] = df['area_ad'].replace(3, 'LI')
df.tail()


# In[125]:


df2 = pd.DataFrame(temp['resp_tg_collapse_ori'], columns=['resp_tg_avg_ori'])
df2['area_tg'] = temp['area_tg']
df2['area_tg'] = df2['area_tg'].replace(1, 'V1')
df2['area_tg'] = df2['area_tg'].replace(2, 'LM')
df2['area_tg'] = df2['area_tg'].replace(3, 'LI')
# df2.tail()


# In[129]:


# df.groupby('area_ad', sort=False).describe().reset_index()


# In[130]:


# df2.groupby('area_tg', sort=False).describe().reset_index()


# In[131]:


# fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharey=False)
# fig.suptitle('resp ad or tg by area')

# # resp ad
# ax = sns.violinplot(ax=axes[0,0], x="area_ad", y="resp_ad", data=df)
# axes[0,0].set_title('resp ad')

# # resp tg
# ax = sns.violinplot(ax=axes[0,1], x="area_tg", y="resp_tg_avg_ori", data=df2)
# axes[0,1].set_title('resp tg')

# # resp ad zoom in
# ax = sns.violinplot(ax=axes[1,0], x="area_ad", y="resp_ad", data=df)
# ax.set(ylim=(-0.1, 0.4))
# axes[1,0].set_title('resp ad zoom in')

# # resp tg zoom in
# ax = sns.violinplot(ax=axes[1,1], x="area_tg", y="resp_tg_avg_ori", data=df2)
# ax.set(ylim=(-0.1, 0.4))
# axes[1,1].set_title('resp tg zoom in')

# fig.tight_layout()


# ### Polar plot of resp_ori by area

# In[134]:


ncell_tg = temp['resp_tg'].shape[0]
nori = temp['resp_tg'].shape[1]
resp_tg_cell = []
resp_tg_sorted = np.zeros((ncell_tg, nori))

for icell in np.arange(ncell_tg):
    resp_tg_cell = temp['resp_tg'][icell, :]
    resp_tg_sorted[icell, :] = np.concatenate(
        (resp_tg_cell[np.argmax(resp_tg_cell):len(resp_tg_cell)], 
         resp_tg_cell[0:np.argmax(resp_tg_cell)])) # put max resp to column zero
    
df3 = pd.DataFrame(resp_tg_sorted)
df3['area'] = temp['area_tg']
df3 = replace_area_name(df3, 'area')
df3.head()


# In[136]:


# df3.groupby('area', sort=False).count().reset_index()


# In[139]:


df_mean = df3.groupby('area', sort=False).mean().reset_index()
df_mean


# In[140]:


df_median = df3.groupby('area', sort=False).median().reset_index()
df_median


# In[141]:


df_std = df3.groupby('area', sort=False).std().reset_index()
df_std


# In[142]:


arr_median = df_median.to_numpy()
arr_median = arr_median[0:3, 1:nori+1]
chop_idx = 4
arr_median = np.hstack((arr_median[:, chop_idx:], arr_median[:, :chop_idx])) # flip to make max resp ori at the middle
arr_median = np.concatenate((arr_median, arr_median[:,0].reshape([3,1])), axis=1)

arr_mean = df_mean.to_numpy()
arr_mean = arr_mean[0:3, 1:nori+1]
arr_mean = np.hstack((arr_mean[:, chop_idx:], arr_mean[:, :chop_idx]))
arr_mean = np.concatenate((arr_mean, arr_mean[:,0].reshape([3,1])), axis=1)

arr_std = df_std.to_numpy()
arr_std = arr_std[0:3, 1:nori+1]
arr_std = np.hstack((arr_std[:, chop_idx:], arr_std[:, :chop_idx]))
arr_std = np.concatenate((arr_std, arr_std[:,0].reshape([3,1])), axis=1)

ncell_area = df3.groupby('area', sort=False).count().reset_index()[0]
ncell_area = ncell_area.to_numpy().reshape([3,1])
arr_sem = arr_std / np.sqrt(ncell_area)

N = nori * 2
theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
radii_V1_med = np.concatenate((arr_median[0], arr_median[0]))
radii_LM_med = np.concatenate((arr_median[1], arr_median[1]))
radii_LI_med = np.concatenate((arr_median[2], arr_median[2]))
radii_V1_mean = np.concatenate((arr_mean[0], arr_mean[0]))
radii_LM_mean = np.concatenate((arr_mean[1], arr_mean[1]))
radii_LI_mean = np.concatenate((arr_mean[2], arr_mean[2]))
width = np.pi / nori


# In[149]:


deg = np.linspace(0.0, 180, 8, endpoint=False)
deg = np.asarray([*deg, 180])
# deg = deg[:-1]

# plt.style.use('seaborn-whitegrid')
# plt.figure(figsize=(6, 5))

# plt.errorbar(deg-1, arr_median[0], yerr=arr_sem[0], label='V1', color='blue', alpha=0.6)
# plt.errorbar(deg, arr_median[1], yerr=arr_sem[1], label='LM', color='orange', alpha=0.6)
# plt.errorbar(deg+1, arr_median[2], yerr=arr_sem[2], label='LI', color='green', alpha=0.6)

# plt.xlim([-15,195])
# plt.ylim([0.02,0.16])
# plt.yticks(np.arange(0.02, 0.18, step=0.02))
# plt.tick_params(
#     axis='x',          # changes apply to the x-axis
#     which='both',      # both major and minor ticks are affected
#     bottom=False,      # ticks along the bottom edge are off
#     top=False,         # ticks along the top edge are off
#     labelbottom=False) # labels along the bottom edge are off
# plt.xlabel('grating orientations aligned')
# plt.ylabel('dF/F')

# figtitle = 'Tuning specificity decreases along ventral stream'
# plt.suptitle(figtitle, size=16)
# plt.legend();
# # plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);

# plt.savefig(root_path + 'plot/poster/' + figtitle + '.pdf', format='pdf')


# In[158]:


# plt.figure(figsize=(15, 10))
# ax = plt.subplot(121, projection='polar')
# ax.bar(theta, radii_V1_med, width=width/1.2, bottom=0.0, color='C0', alpha=0.4)
# ax.bar(theta, radii_LM_med, width=width/1.1, bottom=0.0, color='C1', alpha=0.4)
# ax.bar(theta, radii_LI_med, width=width/1.0, bottom=0.0, color='C3', alpha=0.4)
# plt.yticks(np.arange(0, 0.16, step=0.04))
# ax.legend(['V1 239', 'LM 107', 'LI 66'])
# plt.title('median ori resp')

# ax = plt.subplot(122, projection='polar')
# ax.bar(theta, radii_V1_mean, width=width/1.2, bottom=0.0, color='C0', alpha=0.4)
# ax.bar(theta, radii_LM_mean, width=width/1.1, bottom=0.0, color='C1', alpha=0.4)
# ax.bar(theta, radii_LI_mean, width=width/1.0, bottom=0.0, color='C3', alpha=0.4)
# plt.title('mean ori resp')
# plt.yticks(np.arange(0, 0.18, step=0.04))
# plt.show()


# ## circular variance
# cirvar is calculated with positive values only. all negative responses are rectified to 0.  
# https://stackoverflow.com/questions/52856232/scipy-circular-variance  
# https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Circular_Data_Analysis.pdf  
# https://docs.astropy.org/en/stable/api/astropy.stats.circvar.html

# In[143]:


file_path = os.path.join(root_path, "plot/CV SNR OSI R2 ori_perc by area - why HVA lack well fit", "corr_well_fit_HVA_w_cirvar.mat").replace("\\","/")
temp = loadmat(file_path)

ori_list = np.arange(0, 180, 22.5)
# temp['dfof_tg_ori'].shape

ncell = temp['dfof_tg_ori'].shape[0]
nori = 8
cirvar = np.pi * np.ones((ncell, 1))

for icell in range(0, ncell):
    pseudo = [] # pretend there is a dist of angles, to calculate circular variance
    
    for j in range(0, nori):
        pseudo.append([ori_list[j]] * int(temp['dfof_tg_ori'][icell][j])) 
        pseudo_flat = [item for sublist in pseudo for item in sublist]
        cell_data = np.asarray(pseudo_flat)*u.deg
        cirvar[icell] = float(circvar(cell_data))
        
cirvar.shape


# ## ori_perc correlates with?
# - areacode, vis_driven, well_fit  
# - ori_perc  
# - OSI, R2, width, SNR

# In[144]:


file_path = os.path.join(root_path, "plot/CV SNR OSI R2 ori_perc by area - why HVA lack well fit", "corr_well_fit_HVA.mat").replace("\\","/")
temp = loadmat(file_path)

del temp['__header__']
del temp['__version__']
del temp['__globals__']
temp.keys()


# In[145]:


df = pd.DataFrame(temp['area_merge'], columns=['area'])
for key in temp.keys():
    df[key] = temp[key]

df = df.drop(columns='area_merge')
df['vis_well_fit'] = df['well_fit'] & df['vis'] # true well fit: both visually driven and well fit
df['cirvar'] = cirvar
df['area'] = df['area'].replace(1, 'V1').replace(2, 'LM').replace(3, 'LI')

df = df[['area', 'vis', 'vis_well_fit', 'ori_perc_all', 'well_fit',
    'OSI', 'cirvar', 'sharp', 
    'R2', 'SSE', 'coeff_var']]
df.describe()


# ## corr btw ori_perc vs all

# In[157]:


df['area'] = df['area'].replace('V1',1).replace('LM',2).replace('LI',3)
df_corr = df.corr(method ='pearson')
df_corr[['area', 'vis', 'vis_well_fit', 'ori_perc_all']]


# In[158]:


df_corr


# In[159]:


sns.set_theme(style="white")

f, ax = plt.subplots(figsize=(11, 9))
mask = np.triu(np.ones_like(df_corr, dtype=bool)) # Generate a mask for the upper triangle
cmap = sns.diverging_palette(230, 20, as_cmap=True) # Generate a custom diverging colormap

sns.heatmap(df_corr, mask=mask, cmap=cmap, vmin=-1, vmax=1, annot=True, 
            center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})
plt.xticks(rotation=45);


# In[160]:


df2 = df[['vis_well_fit','ori_perc_all','OSI','R2']]
df_corr_only = df2.corr(method ='pearson')

sns.set_theme(style="white")

f, ax = plt.subplots(figsize=(8,5))
mask = np.triu(np.ones_like(df_corr_only, dtype=bool)) # Generate a mask for the upper triangle
cmap = sns.diverging_palette(230, 20, as_cmap=True) # Generate a custom diverging colormap

sns.heatmap(df_corr_only, mask=mask, cmap=cmap, vmin=-1, vmax=1, annot=True, 
            center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})
plt.xticks(rotation=45); plt.yticks(rotation=45);


# ## find combined metrics of well_fit

# In[161]:


df3 = df[['area','vis','vis_well_fit','ori_perc_all','OSI','R2']]
df3.head()


# In[162]:


plt.figure(figsize=(8,5))
plt.hist(df3[df3.vis_well_fit == 1].R2, bins=20, alpha=0.5, label="well fit");
plt.hist(df3[df3.vis_well_fit == 0].R2, bins=100, alpha=0.5, label="poor fit");
plt.legend();
# plt.xlim(0.5,1);


# In[163]:


plt.figure(figsize=(8,5))
plt.hist(df3[df3.vis_well_fit == 1].OSI, bins=5, alpha=0.7, label="well fit");
plt.hist(df3[df3.vis_well_fit == 0].OSI, bins=100, alpha=0.5, label="poor fit");
plt.legend();
# plt.xlim(0,2);


# In[164]:


df3['re_well_fit'] = ((df3.R2 > 0.8) & (df3.OSI > 0.2) & (df3.OSI < 1.5))
df3['re_vis_well_fit'] = ((df3.R2 > 0.8) & (df3.OSI > 0.2) & (df3.OSI < 1.5) & (df3.vis == 1))
df3.head()


# In[165]:


df3.corr()


# In[166]:


sns.set_theme(style="white")

f, ax = plt.subplots(figsize=(11, 9))
mask = np.triu(np.ones_like(df3.corr(), dtype=bool)) # Generate a mask for the upper triangle
cmap = sns.diverging_palette(230, 20, as_cmap=True) # Generate a custom diverging colormap

sns.heatmap(df3.corr(), mask=mask, cmap=cmap, vmin=-1, vmax=1, annot=True, 
            center=0, square=True, linewidths=.5, cbar_kws={"shrink": .5})
plt.xticks(rotation=45);


# ### overlap btw new vs old well fit cells

# In[167]:


df4 = df3.copy()
df4 = df4.drop(columns=['re_well_fit'])
df4['agree'] = (df3.re_vis_well_fit == df3.vis_well_fit).astype(int)
df4['overlap'] = ((df3.re_vis_well_fit) & (df3.vis_well_fit)).astype(int)
df4.re_vis_well_fit = df4.re_vis_well_fit.astype(int)


# In[168]:


t1 = df4[['area','re_vis_well_fit','vis_well_fit','overlap','agree']].groupby('area', sort=False).mean().reset_index().values
t2 = df4[['area','re_vis_well_fit','vis_well_fit','overlap','agree']].groupby('area', sort=False).count().reset_index().values

df4[['area','re_vis_well_fit','vis_well_fit','overlap','agree']].groupby('area', sort=False).mean().reset_index()


# In[169]:


tt = (t1[:,1:] * t2[:,1:]).astype(int)
df5 = pd.DataFrame(tt)
df5.columns = ['re_vis_well_fit', 'vis_well_fit', 'overlap', 'agree']
df5['total'] = t2[:,1].astype(int)
df5.index = ['V1', 'LM', 'LI']
df5


# ### re_well_fit tuning bias across areas
# filter: vis_driven by adapter & re_well_fit

# In[170]:


file_path = os.path.join(root_path, "mat", "adp_dataset_master.xlsx").replace("\\","/")
meta = pd.read_excel(file_path, index_col=None)

meta = meta[meta.seg == 'segmented']
meta = meta[meta.mouse <= 1324]
meta = meta[['mouse','date','area']]
meta.mouse = meta.mouse.astype(int)
meta.date = meta.date.astype(int)
meta


# In[171]:


nset = len(meta.index); ncell = []; nori = 8; nisi = 3; nframe_trial = 223
dir_name = root_path + 'mat/'

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
    well_fit = np.concatenate((well_fit, cell_prop['well_fit_cell']), axis=0)
    
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

ncell, vis_ad.shape, vis_tg.shape, well_fit.shape, ori_pref.shape, fit_param.shape, dfof_ad.shape, dfof_tg.shape, trace.shape


# In[172]:


df4.head()


# In[173]:


tt = ori_pref.copy()
tt[tt > 90] = np.abs(tt[tt > 90] - 180)
tuning_bias = tt[:,[1,2]] - tt[:,[0]];

ori_pref_bin = tt[:,[0]];
ori_pref_bin[ori_pref_bin < 22.5] = 0; ori_pref_bin[ori_pref_bin > 67.5] = 90; 
ori_pref_bin[(ori_pref_bin >= 22.5) & (ori_pref_bin <= 67.5)] = 45; 


# In[174]:


meta['ncell'] = ncell
mouse_cell = [item for item, count in zip(meta.mouse, meta.ncell) for i in range(count)]
area_cell = [item for item, count in zip(meta.area, meta.ncell) for i in range(count)]
meta_cell = pd.DataFrame({'mouse': mouse_cell, 'area': area_cell})
meta_cell_750 = meta_cell.copy(); meta_cell_750['isi'] = 750
meta_cell_250 = meta_cell.copy(); meta_cell_250['isi'] = 250
meta_cell_isi = pd.concat([meta_cell_750, meta_cell_250], ignore_index=True)

df_adp_tune = meta_cell_isi.copy()
df_adp_tune['tuning_bias'] = tuning_bias.flatten('F')
df_adp_tune['ori_pref_bin'] = np.concatenate((ori_pref_bin, ori_pref_bin), axis=0)

df_adp_tune['vis_ad'] = np.concatenate((vis_ad, vis_ad), axis=0)
df_adp_tune['well_fit'] = np.concatenate((well_fit, well_fit), axis=0)
df_adp_tune['re_well_fit'] = np.concatenate((df4.re_vis_well_fit, df4.re_vis_well_fit), axis=0) 
df_adp_tune = df_adp_tune[ df_adp_tune['vis_ad'] == 1 ]
# df_adp_tune = df_adp_tune[ df_adp_tune['well_fit'] == 1 ]
df_adp_tune = df_adp_tune[ df_adp_tune['re_well_fit'] == True ]

b, c = df_adp_tune.iloc[0].copy(), df_adp_tune.iloc[1].copy() 
df_adp_tune.iloc[0], df_adp_tune.iloc[1] = c, b # swap row 0 & 1 to sort df.gb properly
df_adp_tune.reset_index() #.head(20)


# In[175]:


df_adp_tune[['tuning_bias','area','isi','ori_pref_bin']].groupby(['area'], sort=False).describe()


# In[176]:


df_adp_tune[['tuning_bias','area','isi','ori_pref_bin']].groupby(['area','isi','ori_pref_bin'], sort=False).describe()


# In[177]:


sns.set_style("whitegrid")
plt.figure(figsize=(8,5))
ax = sns.violinplot(data=df_adp_tune[df_adp_tune.area == 'V1'], 
                    x="ori_pref_bin", y="tuning_bias", hue="isi", 
                    split=True, inner="quart", palette="Set3")
sns.despine(left=True)
ax.set(xlabel = '|pref - adapter|', ylabel = 'tuning bias', title = 'V1 tuning bias after adaptation')
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);


# In[178]:


tt = df_adp_tune[['tuning_bias','area','isi','ori_pref_bin']].groupby(['area','isi','ori_pref_bin'], sort=False).mean().reset_index()

ax = sns.lineplot(data=tt, x="ori_pref_bin", y="tuning_bias", hue="area", style="isi", ci='sd');
ax.set(xlabel = '|pref - adapter|', ylabel = 'tuning bias', title = 'mean tuning bias after adaptation')
plt.xlim([-10,100]); plt.xticks([0,45,90])
plt.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.);
plt.hlines(0, 0, 90, linestyles='dotted');


# In[213]:


# bias_mean = df_adp_tune[['tuning_bias','area','isi','ori_pref_bin']].groupby(['area','isi','ori_pref_bin'], sort=False).mean().reset_index().to_numpy()[:,-1]

# bias_sem = df_adp_tune[['tuning_bias','area','isi','ori_pref_bin']].groupby(['area','isi','ori_pref_bin'], sort=False).std().reset_index().to_numpy()[:,-1] / np.sqrt(df_adp_tune[['tuning_bias','area','isi','ori_pref_bin']].groupby(['area','isi','ori_pref_bin'], sort=False).count().reset_index().to_numpy()[:,-1].astype(float))


# ## fig 1C

# In[227]:


sns.set_style("ticks")
sns.set_context("talk", font_scale = 1)
fig, axes = plt.subplots(1, 1, figsize=(6, 5), sharey=True)

ax = sns.violinplot(x="area", y="adp_mag", col="isi",
                data=df_adp_mag_thres[df_adp_mag_thres.isi==250], kind="violin",
                height=4, aspect=1, palette='pastel');

# statistical annotation
x1, x2 = 0, 0.98
y, h, col = np.max(df_adp_mag_thres[df_adp_mag_thres.isi==250].adp_mag.values) + 0.4, 0.3, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, "****", ha='center', va='bottom', color=col, fontsize=16)

x1, x2 = 1.02, 2
y, h, col = np.max(df_adp_mag_thres[df_adp_mag_thres.isi==250].adp_mag.values) + 0.1, 0.3, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, "****", ha='center', va='bottom', color=col, fontsize=16)

x1, x2 = 0, 2
y, h, col = np.max(df_adp_mag_thres[df_adp_mag_thres.isi==250].adp_mag.values) + 1, 0.3, 'k'
plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
plt.text((x1+x2)*.5, y+h, "****", ha='center', va='bottom', color=col, fontsize=16)

plt.xlabel('')
plt.ylabel('adaptation index', fontsize=18)
ax.set_xticklabels(['V1', 'LM', 'LI'], fontsize=18)
ax.set_yticks(np.arange(-3,3.1))
ax.set_ylim((-3,3))

ax.tick_params(axis="x", labelsize=16)
ax.tick_params(axis="y", labelsize=16)
sns.despine(ax=ax)

figtitle = 'Adaptation magnitude increases along ventral stream'
# ax.fig.suptitle(figtitle, size=16)
# ax.fig.subplots_adjust(top=.8)
plt.savefig(root_folder + 'plot/poster/' + figtitle + '.pdf', format='pdf')


# In[228]:


sns.set_theme(style="ticks")
sns.set_context("talk", font_scale = 1)

fig, axes = plt.subplots(1, 1, figsize=(6, 5), sharey=True)
ax = sns.pointplot(ax=axes, x=df_adp_tune.ori_pref_bin.astype(int), y="tuning_bias", hue="area", 
                   data=df_adp_tune[df_adp_tune.isi==250], ci=68, dodge=True,
                   linestyles=["-", "--", "dashdot"])
plt.legend(title='', loc='upper right')
plt.setp(ax.collections, alpha=.5);
plt.setp(ax.lines, alpha=.7);
plt.ylim(-20,20)
plt.xlabel('|preferred orientation - adapter orientation|', labelpad=5, fontsize=18)
plt.ylabel('tuning bias', fontsize=18, labelpad=-8);

ax.tick_params(axis="x", labelsize=16)
ax.tick_params(axis="y", labelsize=16)
sns.despine(ax=ax)

figtitle = 'Tuning bias increases along ventral stream'
# fig.suptitle(figtitle)
plt.tight_layout()
plt.savefig(root_path + 'plot/poster/' + figtitle + '.pdf', format='pdf')


# In[229]:


# sns.set_style("ticks")
# sns.set_context("talk", font_scale = 1)

# fig, ax = plt.subplots(1, 2)

# # ax[0].plot(range(10), 'r') #row=0, col=0


# ax = sns.violinplot(x="area", y="adp_mag", col="isi",
#                 data=df_adp_mag_thres[df_adp_mag_thres.isi==250], kind="violin",
#                 height=4, aspect=1, palette='pastel', ax=ax[0]);

# # statistical annotation
# x1, x2 = 0, 0.98
# y, h, col = np.max(df_adp_mag_thres[df_adp_mag_thres.isi==250].adp_mag.values) + 0.4, 0.3, 'k'
# plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
# plt.text((x1+x2)*.5, y+h, "****", ha='center', va='bottom', color=col, fontsize=16)

# x1, x2 = 1.02, 2
# y, h, col = np.max(df_adp_mag_thres[df_adp_mag_thres.isi==250].adp_mag.values) + 0.1, 0.3, 'k'
# plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
# plt.text((x1+x2)*.5, y+h, "****", ha='center', va='bottom', color=col, fontsize=16)

# x1, x2 = 0, 2
# y, h, col = np.max(df_adp_mag_thres[df_adp_mag_thres.isi==250].adp_mag.values) + 1, 0.3, 'k'
# plt.plot([x1, x1, x2, x2], [y, y+h, y+h, y], lw=1.5, c=col)
# plt.text((x1+x2)*.5, y+h, "****", ha='center', va='bottom', color=col, fontsize=16)

# plt.xlabel('')
# plt.ylabel('adaptation index', fontsize=18)
# ax.set_xticklabels(['V1', 'LM', 'LI'], fontsize=18)
# ax.set_yticks(np.arange(-3,3.1))
# ax.set_ylim((-3,3))

# ax.tick_params(axis="x", labelsize=16)
# ax.tick_params(axis="y", labelsize=16)
# sns.despine(ax=ax)


# # ax[1].plot(range(10), 'b') #row=1, col=0
# plt.show()


# In[ ]:




