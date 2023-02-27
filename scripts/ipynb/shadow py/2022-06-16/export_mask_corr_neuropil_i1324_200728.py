#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import scipy.io
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
# https://github.com/MouseLand/suite2p/issues/166
# https://suite2p.readthedocs.io/_/downloads/en/latest/pdf/


# In[11]:


STATS_FILENAME = 'D:\Lan_temp\suite2p\plane0_i1324_200728_003/stat.npy'
st = np.load(STATS_FILENAME, allow_pickle=True)
ops_FILENAME = 'D:\Lan_temp\suite2p\plane0_i1324_200728_003/ops.npy'
ops = np.load(ops_FILENAME, allow_pickle=True).item()

n_rois = st.shape[0]
height = 264; width = 796 # Image size in pixels
mask = np.zeros((n_rois, height, width))

for i in range(n_rois):
    mask[i, st[i]['ypix'], st[i]['xpix']] = st[i]['lam']
'lam'.upper()


# In[4]:


iscell = np.load('D:\Lan_temp\suite2p\plane0_i1324_200728_003/iscell.npy', allow_pickle=True)[:, 0].astype(bool)
iscell.shape, mask.shape # find out why mask is not digital


# In[6]:


stat = st
im = np.zeros((ops['Ly'], ops['Lx']))

for n in range(0,st.shape[0]):
    ypix = stat[n]['ypix'][~stat[n]['overlap']]
    xpix = stat[n]['xpix'][~stat[n]['overlap']]
    im[ypix,xpix] = n+1
    
plt.figure(figsize=(20,10));
plt.imshow(im)
plt.show()


# In[7]:


plt.figure(figsize=(20,10));
plt.imshow(mask[100])


# In[8]:


mask_all = sum(mask,1)
plt.figure(figsize=(20,10))
plt.imshow(mask_all)


# In[10]:


plt.hist(mask[100]);
plt.xlim([0,0.004]);


# In[9]:


sum(sum(mask[100] > 0.0035))


# In[7]:


i=1
mask[i, st[i]['ypix'], st[i]['xpix']] 


# In[9]:


st[i]['ypix']


# In[9]:


st[i]['xpix']


# # corr matrix

# In[5]:


cd //duhs-user-nc1.dhe.duke.edu/dusom_glickfeldlab/All_Staff/home/lan/Analysis/2P/200728_i1324/200728_i1324_runs-003


# In[8]:


ref = scipy.io.loadmat('200728_i1324_runs-003_TCs_addfake')
ref_npSub_tc = np.asarray(ref['npSub_tc'])
ref_np_tc = np.asarray(ref['np_tc'])

ref1 = pd.DataFrame(ref_np_tc)
ref2 = pd.DataFrame(ref_npSub_tc)
corr1 = ref1.corr()
corr2 = ref2.corr() # after np sub

ref_npSub_tc.shape, ref_np_tc.shape


# In[56]:


fig, axes = plt.subplots(1, 2, figsize=(20, 5), sharey=True);

sns.heatmap(corr1, ax=axes[0], vmin=-0.1, vmax=1.0)
sns.heatmap(corr2, ax=axes[1], vmin=-0.1, vmax=1.0)

plt.show()


# In[65]:


fig, axes = plt.subplots(1, 2, figsize=(15, 5), sharex=True, sharey=True);
nbin = 20

axes[0].hist(corr1.to_numpy().flatten(), bins=nbin)
axes[0].axvline(x=np.median(corr1.to_numpy().flatten()), linewidth=2, color='g')
axes[0].axvline(x=np.mean(corr1.to_numpy().flatten()), linewidth=2, color='r')
axes[1].hist(corr2.to_numpy().flatten(), bins=nbin)
axes[1].axvline(x=np.median(corr2.to_numpy().flatten()), linewidth=2, color='g')
axes[1].axvline(x=np.mean(corr2.to_numpy().flatten()), linewidth=2, color='r')

plt.tight_layout()
plt.show()


# In[58]:


f_original = np.load('D:\Lan_temp\suite2p\plane0_i1324_200728_003\F.npy', allow_pickle=True)
f_neuropils = np.load('D:\Lan_temp\suite2p\plane0_i1324_200728_003\Fneu.npy', allow_pickle=True)

ops['neucoeff'] = 0.7 # default setting
f_neuron = f_original - ops['neucoeff'] * f_neuropils

f_neuron = f_neuron[iscell,:]
f_ori = f_original[iscell,:]
f_ori.shape, f_neuron.shape


# In[32]:


s1 = pd.DataFrame(f_ori.transpose())
s2 = pd.DataFrame(f_neuron.transpose())
corr1s = s1.corr()
corr2s = s2.corr() # after np sub

fig = plt.figure(figsize=[20,5])
plt.subplot(1, 2, 1)
sns.heatmap(corr1s)

plt.subplot(1, 2, 2)
sns.heatmap(corr2s)
plt.show()


# In[39]:


fig, axes = plt.subplots(1, 2, figsize=(20, 5), sharey=True);

sns.heatmap(corr1s, ax=axes[0], vmin=-0.45, vmax=1.0)
sns.heatmap(corr2s, ax=axes[1], vmin=-0.45, vmax=1.0)

plt.show()


# In[75]:


fig, axes = plt.subplots(1, 2, figsize=(15, 5), sharex=True, sharey=True);
nbin = 20

axes[0].hist(corr1s.to_numpy().flatten(), bins=nbin)
axes[0].axvline(x=np.nanmedian(corr1s.to_numpy().flatten()), linewidth=2, color='g')
axes[0].axvline(x=np.nanmean(corr1s.to_numpy().flatten()), linewidth=2, color='r')
axes[1].hist(corr2s.to_numpy().flatten(), bins=nbin)
axes[1].axvline(x=np.median(corr2s.to_numpy().flatten()), linewidth=2, color='g')
axes[1].axvline(x=np.mean(corr2s.to_numpy().flatten()), linewidth=2, color='r')

plt.tight_layout()
plt.show()


# In[79]:


ops['max_iterations']

