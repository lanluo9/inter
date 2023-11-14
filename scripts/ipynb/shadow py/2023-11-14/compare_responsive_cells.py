#!/usr/bin/env python
# coding: utf-8

# In[48]:


from pathlib import Path
from tempfile import TemporaryDirectory

import matplotlib.pyplot as plt
import numpy as np
import scipy.io

import suite2p
from suite2p import run_s2p, default_ops
# ops = default_ops() # populates ops with the default options
import sys
sys.path.insert(0, 'C:/Users/lan/Documents/repos/suite2p') # option to import from github folder

import seaborn as sns
import matplotlib as mpl
mpl.rcParams.update({
    'axes.spines.left': True,
    'axes.spines.bottom': False,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'legend.frameon': False,
    'figure.subplot.wspace': .01,
    'figure.subplot.hspace': .01,
#     'figure.figsize': (18, 13),
    'ytick.major.left': True,
})
jet = mpl.cm.get_cmap('jet')
jet.set_bad(color='k')


# In[2]:


t = 'D:\Lan_temp\suite2p\plane0_i1322_200803_002'
output_op = np.load(Path(t).joinpath('ops.npy'), allow_pickle=True).item()
output_op['save_path'] = t


# In[6]:


stats_file = Path(output_op['save_path']).joinpath('stat.npy')
stats = np.load(stats_file, allow_pickle=True)
iscell = np.load(Path(output_op['save_path']).joinpath('iscell.npy'), allow_pickle=True)[:, 0].astype(bool)
iscell.shape


# In[81]:


f_original = np.load(Path(output_op['save_path']).joinpath('F.npy'))
f_neuropils = np.load(Path(output_op['save_path']).joinpath('Fneu.npy'))
spks = np.load(Path(output_op['save_path']).joinpath('spks.npy'))
f_cells = f_original - ops['neucoeff'] * f_neuropils


# In[86]:


# plt.figure(figsize=[20,20])
rois = np.arange(len(f_cells))[::2000]
frame_range = np.arange(100,400)

for i, roi in enumerate(rois):
    plt.subplot(len(rois), 1, i+1, )
    f = f_cells[roi][frame_range]
    f_neu = f_neuropils[roi][frame_range]
    sp = spks[roi][frame_range]
    
    # Adjust spks range to match range of fluroescence traces
    fmax = np.maximum(f.max(), f_neu.max())
    fmin = np.minimum(f.min(), f_neu.min())
    frange = fmax - fmin 
    sp /= sp.max()
    sp *= frange
    
    plt.plot(f, label="Cell Fluorescence")
    plt.plot(f_neu, label="Neuropil Fluorescence")
    plt.plot(sp + fmin, label="Deconvolved")
#     plt.xticks(np.arange(0, f_cells.shape[1], f_cells.shape[1]/10))
#     plt.ylabel(f"ROI {roi}", rotation=0)
    plt.xlabel("frame")
    if i == 0:
        plt.legend(bbox_to_anchor=(0.93, 2))


# In[25]:


sum(iscell), sum(~iscell)


# In[50]:


f_neuron = f_cells[iscell,:]
# f_neuron = f_original[iscell,:]
f_neuron.shape


# In[51]:


tt = np.mean(f_neuron, axis=0)
tt.max(), tt.min(), np.median(tt)


# In[52]:


sns.histplot(tt)


# In[35]:


plt.figure(figsize=[15,5])
plt.plot(tt[0:1000])
# plt.yticks(np.arange(0, 1, step=0.2))


# In[53]:


temp_dic = {"f_neuron": f_neuron}
scipy.io.savemat("f_neuron.mat", temp_dic)


# ### after selecting responsive cells in matlab

# In[57]:


stats


# In[72]:


# stat = np.load('stat.npy')
# ops = np.load('ops.npy').item()

stat = stats
ops = output_op
im = np.zeros((ops['Ly'], ops['Lx']))

ncells = len(iscell)
for n in range(0,ncells):
    ypix = stat[n]['ypix'][~stat[n]['overlap']]
    xpix = stat[n]['xpix'][~stat[n]['overlap']]
    im[ypix,xpix] = n+1

plt.imshow(im)
plt.show()


# In[78]:


stat[0]['ypix'][~stat[0]['overlap']].shape

