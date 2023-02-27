#!/usr/bin/env python
# coding: utf-8

# In[122]:


import pandas as pd
import numpy as np
import scipy
from scipy.io import loadmat
from sklearn.linear_model import LinearRegression
import seaborn as sns
import matplotlib.pyplot as plt
from astropy.stats import circvar
from astropy import units as u


# In[123]:


def replace_area_name(df, area_key):
# df = dataframe containing an area array column
# area_key = string

    df[area_key] = df[area_key].replace(1, 'V1')
    df[area_key] = df[area_key].replace(2, 'LM')
    df[area_key] = df[area_key].replace(3, 'LI')
    
    return df


# ### OSI by area
# #### for only well-fit cells

# In[27]:


OSI_area = loadmat('C:/Users/lan/Documents/repos/inter/plot/OSI_area.mat')
df = pd.DataFrame(OSI_area['OSI'], columns=['OSI'])
df['area'] = OSI_area['area_merge']
df['area'] = df['area'].replace(1, 'V1')
df['area'] = df['area'].replace(2, 'LM')
df['area'] = df['area'].replace(3, 'LI')
df.tail()


# In[33]:


df.groupby('area', sort=False).describe().reset_index()


# In[30]:


ax = sns.violinplot(x="area", y="OSI", data=df)


# #### OSI by area for all vis driven (by tg) cells

# In[246]:


OSI_area = loadmat('C:/Users/lan/Documents/repos/inter/plot/OSI_area_vis.mat')
df = pd.DataFrame(OSI_area['OSI_all'], columns=['OSI'])
df['well_fit'] = OSI_area['well_fit']
df['area'] = OSI_area['area']
df = replace_area_name(df, 'area')
df.tail()


# In[247]:


df.groupby('well_fit', sort=False).describe().reset_index()


# In[257]:


df[['area','well_fit']].groupby('area', sort=False).count().reset_index()


# In[250]:


# plt.figure(figsize=(15, 10))
ax = sns.violinplot(x="area", y="OSI", hue="well_fit",
                    data=df, palette="Set3", split=True)


# ### Circular Variance by area
# 
# https://stackoverflow.com/questions/52856232/scipy-circular-variance  
# https://ncss-wpengine.netdna-ssl.com/wp-content/themes/ncss/pdf/Procedures/NCSS/Circular_Data_Analysis.pdf  
# https://docs.astropy.org/en/stable/api/astropy.stats.circvar.html

# In[48]:


temp = loadmat('C:/Users/lan/Documents/repos/inter/plot/CirVar_area.mat')
ori_list = np.arange(0, 180, 22.5)
temp['dfof_tg_ori']


# In[170]:


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


# In[172]:


cirvar.shape


# In[174]:


df2 = pd.DataFrame(cirvar, columns=['cirvar'])
df2['area'] = temp['area_merge']
df2['area'] = df['area'].replace(1, 'V1')
df2['area'] = df['area'].replace(2, 'LM')
df2['area'] = df['area'].replace(3, 'LI')
df2.tail()


# In[175]:


df2.groupby('area', sort=False).describe().reset_index()


# In[177]:


ax = sns.violinplot(x="area", y="cirvar", data=df2)


# In[184]:


df3 = pd.DataFrame(df2['area'], columns=['area'])
df3['cir_converge'] = 1 - df2['cirvar']
df3.groupby('area', sort=False).describe().reset_index()


# In[186]:


ax = sns.violinplot(x="area", y="cir_converge", data=df3)


# ### Coefficient of variance (reverse SNR) by area

# In[188]:


temp = loadmat('C:/Users/lan/Documents/repos/inter/plot/CV_area_well_fit.mat')
df4 = pd.DataFrame(temp['coeff_var'], columns=['coeff_var'])
df4['area'] = temp['area_merge']
df4['area'] = df4['area'].replace(1, 'V1')
df4['area'] = df4['area'].replace(2, 'LM')
df4['area'] = df4['area'].replace(3, 'LI')
df4.tail()


# In[189]:


df4.groupby('area', sort=False).describe().reset_index()


# In[192]:


ax = sns.violinplot(x="area", y="coeff_var", data=df4)


# In[2]:


temp = loadmat('C:/Users/lan/Documents/repos/inter/plot/CV_area.mat')
df5 = pd.DataFrame(temp['CV'], columns=['CV'])
df5['well_fit'] = temp['well_fit']
df5['area'] = temp['area']

df5['area'] = df5['area'].replace(1, 'V1')
df5['area'] = df5['area'].replace(2, 'LM')
df5['area'] = df5['area'].replace(3, 'LI')
df5.head()


# In[3]:


df5.groupby('area', sort=False).describe().reset_index()


# In[9]:


df5.groupby(['well_fit'], sort=False).describe().reset_index()


# In[4]:


# plt.figure(figsize=(15, 10))
ax = sns.violinplot(x="area", y="CV", hue="well_fit",
                    data=df5, palette="Set3", split=True)


# In[5]:


plt.figure(figsize=(8, 4))
ax = sns.violinplot(x="area", y="CV", hue="well_fit",
                    data=df5, palette="Set3", split=True)
ax.set(ylim=(-1, 5))


# ### Response Amplitude by area

# In[6]:


temp = loadmat('C:/Users/lan/Documents/repos/inter/plot/amp_area.mat')
df = pd.DataFrame(temp['resp_ad'], columns=['resp_ad'])
df['area_ad'] = temp['area_ad']
df['area_ad'] = df['area_ad'].replace(1, 'V1')
df['area_ad'] = df['area_ad'].replace(2, 'LM')
df['area_ad'] = df['area_ad'].replace(3, 'LI')
df.tail()


# In[12]:


df2 = pd.DataFrame(temp['resp_tg_collapse_ori'], columns=['resp_tg_avg_ori'])
df2['area_tg'] = temp['area_tg']
df2['area_tg'] = df2['area_tg'].replace(1, 'V1')
df2['area_tg'] = df2['area_tg'].replace(2, 'LM')
df2['area_tg'] = df2['area_tg'].replace(3, 'LI')
df2.tail()


# In[13]:


df.groupby('area_ad', sort=False).describe().reset_index()


# In[14]:


df2.groupby('area_tg', sort=False).describe().reset_index()


# In[28]:


fig, axes = plt.subplots(2, 2, figsize=(15, 10), sharey=False)
fig.suptitle('resp ad or tg by area')

# resp ad
ax = sns.violinplot(ax=axes[0,0], x="area_ad", y="resp_ad", data=df)
# axes[0,0].set_title('resp ad')

# resp tg
ax = sns.violinplot(ax=axes[0,1], x="area_tg", y="resp_tg_avg_ori", data=df2)
# axes[0,1].set_title('resp tg')

# resp ad zoom in
ax = sns.violinplot(ax=axes[1,0], x="area_ad", y="resp_ad", data=df)
ax.set(ylim=(-0.1, 0.4))
# axes[1,0].set_title('resp ad zoom in')

# resp tg zoom in
ax = sns.violinplot(ax=axes[1,1], x="area_tg", y="resp_tg_avg_ori", data=df2)
ax.set(ylim=(-0.1, 0.4))
# axes[1,1].set_title('resp tg zoom in')


# ### Polar plot of resp_ori by area

# In[79]:


ncell_tg = temp['resp_tg'].shape[0]
nori = temp['resp_tg'].shape[1]
resp_tg_cell = []
resp_tg_sorted = np.zeros((ncell_tg, nori))


# In[117]:


for icell in np.arange(ncell_tg):
    resp_tg_cell = temp['resp_tg'][icell, :]
    resp_tg_sorted[icell, :] = np.concatenate((resp_tg_cell[np.argmax(resp_tg_cell):len(resp_tg_cell)], resp_tg_cell[0:np.argmax(resp_tg_cell)]))


# In[132]:


df3 = pd.DataFrame(resp_tg_sorted)
df3['area'] = temp['area_tg']
df3 = replace_area_name(df3, 'area')
df3


# In[135]:


df3.groupby('area', sort=False).count().reset_index()


# In[139]:


df_median = df3.groupby('area', sort=False).median().reset_index()
df_median


# In[140]:


df_mean = df3.groupby('area', sort=False).mean().reset_index()
df_mean


# In[138]:


df3.groupby('area', sort=False).std().reset_index()


# In[214]:


arr_median = df_median.to_numpy()
arr_median = arr_median[0:3, 1:nori+1]

arr_mean = df_mean.to_numpy()
arr_mean = arr_mean[0:3, 1:nori+1]

N = nori * 2
theta = np.linspace(0.0, 2 * np.pi, N, endpoint=False)
radii_V1_med = np.concatenate((arr_median[0], arr_median[0]))
radii_LM_med = np.concatenate((arr_median[1], arr_median[1]))
radii_LI_med = np.concatenate((arr_median[2], arr_median[2]))
radii_V1_mean = np.concatenate((arr_mean[0], arr_mean[0]))
radii_LM_mean = np.concatenate((arr_mean[1], arr_mean[1]))
radii_LI_mean = np.concatenate((arr_mean[2], arr_mean[2]))
width = np.pi / nori


# In[233]:


plt.figure(figsize=(15, 10))
ax = plt.subplot(121, projection='polar')
ax.bar(theta, radii_V1_med, width=width/1.2, bottom=0.0, color='C0', alpha=0.4)
ax.bar(theta, radii_LM_med, width=width/1.1, bottom=0.0, color='C1', alpha=0.4)
ax.bar(theta, radii_LI_med, width=width/1.0, bottom=0.0, color='C3', alpha=0.4)
plt.yticks(np.arange(0, 0.16, step=0.04))
ax.legend(['V1 239', 'LM 107', 'LI 66'])

ax = plt.subplot(122, projection='polar')
ax.bar(theta, radii_V1_mean, width=width/1.2, bottom=0.0, color='C0', alpha=0.4)
ax.bar(theta, radii_LM_mean, width=width/1.1, bottom=0.0, color='C1', alpha=0.4)
ax.bar(theta, radii_LI_mean, width=width/1.0, bottom=0.0, color='C3', alpha=0.4)
plt.yticks(np.arange(0, 0.18, step=0.04))
plt.show()


# ### R2/SSE by area

# In[5]:


temp = loadmat('C:/Users/lan/Documents/repos/inter/plot/R2_SSE_area.mat')
df = pd.DataFrame(temp['area'], columns=['area'])
df['R2'] = temp['R2']
df['SSE'] = temp['SSE']
replace_area_name(df, 'area')


# In[10]:


df[['area','R2']].groupby('area', sort=False).describe().reset_index()


# In[11]:


df[['area','SSE']].groupby('area', sort=False).describe().reset_index()


# In[7]:


ax = sns.violinplot(x="area", y="R2", data=df)


# In[8]:


ax = sns.violinplot(x="area", y="SSE", data=df)


# In[117]:


temp = loadmat('C:/Users/lan/Documents/repos/inter/plot/R2_SSE_all_area.mat')
df = pd.DataFrame(temp['area'], columns=['area'])
df['R2'] = temp['R2']
df['SSE'] = temp['SSE']
replace_area_name(df, 'area')


# In[118]:


df[['area','R2']].groupby('area', sort=False).describe().reset_index()


# In[119]:


df[['area','SSE']].groupby('area', sort=False).describe().reset_index()


# In[120]:


ax = sns.violinplot(x="area", y="R2", data=df)


# In[121]:


ax = sns.violinplot(x="area", y="SSE", data=df)


# ## Stableness of fit (ori_perc) by area

# In[14]:


temp = loadmat('C:/Users/lan/Documents/repos/inter/plot/ori_perc_area.mat')
df = pd.DataFrame(temp['area'], columns=['area'])
df['ori_perc'] = temp['ori_perc']
replace_area_name(df, 'area')

OSI_area = loadmat('C:/Users/lan/Documents/repos/inter/plot/OSI_area.mat')
df['OSI'] = OSI_area['OSI']
df.tail()


# In[15]:


df2 = df[['area', 'ori_perc']]
df2.tail()


# In[16]:


df2.groupby('area', sort=False).describe().reset_index()


# ### CDF of ori_perc dist of all cells across areas

# In[62]:


temp = loadmat('C:/Users/lan/Documents/repos/inter/plot/ori_perc_all_area.mat')
df = pd.DataFrame(temp['area'], columns=['area'])
df['ori_perc_all'] = temp['ori_perc_all']
replace_area_name(df, 'area')
df.tail()


# In[23]:


df.groupby('area', sort=False).describe().reset_index()


# In[37]:


# df_perc = df.groupby('area', sort=False).quantile(np.linspace(0,1,101)).copy()
# df_perc.tail()


# In[72]:


data1 = df['ori_perc_all'][df['area'] == 'V1']; data_sorted1 = np.sort(data1)
data2 = df['ori_perc_all'][df['area'] == 'LM']; data_sorted2 = np.sort(data2)
data3 = df['ori_perc_all'][df['area'] == 'LI']; data_sorted3 = np.sort(data3)

p1 = 1. * np.arange(len(data1)) / (len(data1) - 1) # calculate the proportional values of samples
p2 = 1. * np.arange(len(data2)) / (len(data2) - 1)
p3 = 1. * np.arange(len(data3)) / (len(data3) - 1)

fig, axes = plt.subplots(1, 3, figsize=(10, 5), sharey=True)
fig.suptitle('CDF of ori_perc')

ax = axes[0].plot(data_sorted1, p1)
axes[0].set_title('V1')
ax = axes[1].plot(data_sorted2, p2)
axes[1].set_title('LM')
ax = axes[2].plot(data_sorted3, p3)
axes[2].set_title('LI')

# fig = plt.figure()
# ax1 = fig.add_subplot(131)
# ax1.plot(data_sorted1, p1)
# ax1.set_xlabel('$x$')
# ax1.set_ylabel('$p$')

# ax2 = fig.add_subplot(132)
# ax2.plot(data_sorted2, p2)
# ax3 = fig.add_subplot(133)
# ax3.plot(data_sorted3, p3)


# In[78]:


fig = plt.figure()
plt.plot(data_sorted1, p1)
plt.plot(data_sorted2, p2)
plt.plot(data_sorted3, p3)
plt.legend(['V1','LM','LI'])
plt.xlabel("ori_perc")
plt.ylabel("prob")


# ### ori_perc vs OSI

# In[124]:


# for well_fit cells
temp = loadmat('C:/Users/lan/Documents/repos/inter/plot/ori_perc_area.mat')
df = pd.DataFrame(temp['area'], columns=['area'])
df['ori_perc'] = temp['ori_perc']
# replace_area_name(df, 'area')

OSI_area = loadmat('C:/Users/lan/Documents/repos/inter/plot/OSI_area.mat')
df['OSI'] = OSI_area['OSI']
df.tail()


# In[91]:


# rng = np.random.RandomState(0)
x = df['ori_perc']
y = df['OSI']
colors = 4 - df['area']
sizes = 100

plt.figure(figsize=(10,8))
plt.scatter(x, y, c=colors, s=sizes, alpha=0.5, cmap='viridis')
# plt.colorbar();
plt.xlabel("ori_perc of well-fit cells")
plt.ylabel("OSI")


# In[ ]:


# for all vis cells


# In[89]:


X = df['ori_perc'].copy().to_numpy().reshape(-1, 1)
len(X)


# In[101]:


X = df['ori_perc'][df['area'] == 1].copy().to_numpy().reshape(-1, 1)
y = df['OSI'][df['area'] == 1].copy().to_numpy().reshape(-1, 1)
reg = LinearRegression().fit(X, y)
reg.score(X, y)


# In[102]:


X = df['ori_perc'].copy().to_numpy().reshape(-1, 1)
y = df['OSI'].copy().to_numpy().reshape(-1, 1)
reg = LinearRegression().fit(X, y)
reg.score(X, y)


# In[112]:


X = df['ori_perc'][df['area'] == 1].copy().to_numpy()
y = df['OSI'][df['area'] == 1].copy().to_numpy()
r, p = scipy.stats.pearsonr(X, y)
r


# In[113]:


p


# In[114]:


np.corrcoef(X, y)


# In[125]:


X = df['ori_perc'].copy().to_numpy()
y = df['OSI'].copy().to_numpy()
r, p = scipy.stats.pearsonr(X, y)
r


# In[126]:


p


# In[ ]:




