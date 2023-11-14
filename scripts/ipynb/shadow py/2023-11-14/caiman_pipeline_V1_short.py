#!/usr/bin/env python
# coding: utf-8

# ## Prep

# In[1]:


import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import bokeh.plotting as bpl
import cv2
import glob
import logging
import matplotlib.pyplot as plt
plt.style.use('seaborn-talk')
import numpy as np
import os

try:
    cv2.setNumThreads(0)
except():
    pass

try:
    if __IPYTHON__:
        # this is used for debugging purposes only. allows to reload classes
        # when changed
        get_ipython().magic('load_ext autoreload')
        get_ipython().magic('autoreload 2')
except NameError:
    pass

import caiman as cm
from caiman.motion_correction import MotionCorrect
from caiman.source_extraction.cnmf import cnmf as cnmf
from caiman.source_extraction.cnmf import params as params
from caiman.utils.utils import download_demo
from caiman.utils.visualization import plot_contours, nb_view_patches, nb_plot_contour
bpl.output_notebook()

logging.basicConfig(format=
                    "%(relativeCreated)12d [%(filename)s:%(funcName)20s():%(lineno)s] [%(process)d] %(message)s",
#                     filename="/tmp/caiman.log",
                    level=logging.INFO)


# In[2]:


caiman_data_path = 'C:/Users/ll357/Documents/CaImAn/demos/temp_data/'
tif_files = glob.glob(caiman_data_path + 'i44415_220318/' + "*.tif", recursive = True)
print(tif_files)

fnames = tif_files[1]
# fnames = tif_files[-1]
fnames


# In[3]:


param_area = 'V1'
param_ca = 'gcamp8f'
# param_stim = 'grating'
param_stim = 'bunnytop high res'

fr = 30                             # imaging rate in frames per second
decay_time = 0.3                    # length of a typical transient in seconds
# decay_time is an approximation of the time scale over which to expect a significant shift in the calcium signal during a transient. 
# It defaults to 0.4, which is appropriate for fast indicators (GCaMP6f), slow indicators might use 1 or even more. 
# However, decay_time does not have to precisely fit the data, approximations are enough.

# decay_time = 0.8 for gcamp6s, ~750 ms decay time in V1 L2/3: https://elifesciences.org/articles/51675
# V1 layer 2/3: rise and decay time constants for GCaMP6f are 50-100 ms and 200-300 ms; 
# and for GCaMP6s 150-200 ms and ~750 ms (4-5 APs, Figure 3C). https://elifesciences.org/articles/51675

# decay_time = 0.3 for gcamp8f, https://www.janelia.org/jgcamp8-calcium-indicators
# Half-decay time (ms) = 67.4±11.2, should be faster than 6f whose decay_time = 0.4

# motion correction parameters
strides = (48, 48)          # start a new patch for pw-rigid motion correction every x pixels
overlaps = (24, 24)         # overlap between pathes (size of patch strides+overlaps)
max_shifts = (6, 6)          # maximum allowed rigid shifts (in pixels)
max_deviation_rigid = 3     # maximum shifts deviation allowed for patch with respect to rigid shifts
pw_rigid = True             # flag for performing non-rigid motion correction

# parameters for source extraction and deconvolution
p = 1                       # order of the autoregressive system
gnb = 2                     # number of global background components
merge_thr = 0.80            # merging threshold, max correlation allowed # lowered value to reduce fake neurons
rf = 15                     # half-size of the patches in pixels. e.g., if rf=25, patches are 50x50
stride_cnmf = 6             # amount of overlap between the patches in pixels
K = 4                       # number of components per patch
gSig = [8, 8]               # expected half size of neurons in pixels
# gSig kept 4 according to snipaste pixel measurement 
# -> changed to 6 to match with mag 2.8x, and changed K from 4 to 2. 2022-03-08
# -> changed K from 2 to 4 in case it meant splitting patch smaller? 2022-03-26
# -> changed to 8 since gSig=6 was still segmented too small. 8 looks good.
# adjust cell size bc of zoom change around 2021-09

method_init = 'greedy_roi'  # initialization method (if analyzing dendritic data using 'sparse_nmf')
ssub = 1                    # spatial subsampling during initialization
tsub = 1                    # temporal subsampling during intialization

# parameters for component evaluation
min_SNR = 2.0               # signal to noise ratio for accepting a component
rval_thr = 0.90             # space correlation threshold for accepting a component # adjusted
rval_lowest = 0.1           # added to eliminate fake overlapping neurons, default -1 # adjusted higher
cnn_thr = 0.99              # threshold for CNN based classifier
cnn_lowest = 0.50           # neurons with cnn probability lower than this value are rejected # adjusted


opts_dict = {'fnames': fnames,
            'fr': fr,
            'decay_time': decay_time,
            'strides': strides,
            'overlaps': overlaps,
            'max_shifts': max_shifts,
            'max_deviation_rigid': max_deviation_rigid,
            'pw_rigid': pw_rigid,
            'p': p,
            'nb': gnb,
            'rf': rf,
            'K': K, 
            'stride': stride_cnmf,
            'method_init': method_init,
            'rolling_sum': True,
            'only_init': True,
            'ssub': ssub,
            'tsub': tsub,
            'merge_thr': merge_thr, 
            'min_SNR': min_SNR,
            'rval_thr': rval_thr,
            'use_cnn': True,
            'min_cnn_thr': cnn_thr,
            'cnn_lowest': cnn_lowest}

opts = params.CNMFParams(params_dict=opts_dict)


# In[4]:


# param_area = 'LM'
# param_ca = 'gcamp6s'
# # param_stim = 'grating'
# param_stim = 'bunnytop high res slow'

# fr = 30                             # imaging rate in frames per second
# decay_time = 0.8                    # length of a typical transient in seconds
# # revised to account for gcamp6s, ~750 ms decay time in V1 L2/3: https://elifesciences.org/articles/51675
# # decay_time: Length of a typical transient in seconds. decay_time is an approximation of the time scale over which to expect a significant shift in the calcium signal during a transient. It defaults to 0.4, which is appropriate for fast indicators (GCaMP6f), slow indicators might use 1 or even more. However, decay_time does not have to precisely fit the data, approximations are enough.
# # V1 layer 2/3: rise and decay time constants for GCaMP6f are 50-100 ms and 200-300 ms; and for GCaMP6s 150-200 ms and ~750 ms (4-5 APs, Figure 3C). https://elifesciences.org/articles/51675

# # motion correction parameters
# strides = (48, 48)          # start a new patch for pw-rigid motion correction every x pixels
# overlaps = (24, 24)         # overlap between pathes (size of patch strides+overlaps)
# max_shifts = (6,6)          # maximum allowed rigid shifts (in pixels)
# max_deviation_rigid = 3     # maximum shifts deviation allowed for patch with respect to rigid shifts
# pw_rigid = True             # flag for performing non-rigid motion correction

# # parameters for source extraction and deconvolution
# p = 1                       # order of the autoregressive system
# gnb = 2                     # number of global background components
# merge_thr = 0.85            # merging threshold, max correlation allowed
# rf = 15                     # half-size of the patches in pixels. e.g., if rf=25, patches are 50x50
# stride_cnmf = 6             # amount of overlap between the patches in pixels
# K = 2                       # number of components per patch
# gSig = [8, 8]               # expected half size of neurons in pixels
# # gSig kept 4 according to snipaste pixel measurement 
# # -> changed to 6 to match with mag 2.8x, and changed K from 4 to 2. 2022-03-08
# # -> changed to 8 since gSig=6 was still segmented too small. 8 looks good.
# # adjust cell size bc of zoom change around 2021-09

# method_init = 'greedy_roi'  # initialization method (if analyzing dendritic data using 'sparse_nmf')
# ssub = 1                    # spatial subsampling during initialization
# tsub = 1                    # temporal subsampling during intialization

# # parameters for component evaluation
# min_SNR = 2.0               # signal to noise ratio for accepting a component
# rval_thr = 0.90             # space correlation threshold for accepting a component # adjusted
# rval_lowest = 0.0           # added to eliminate fake overlapping neurons, default -1
# cnn_thr = 0.99              # threshold for CNN based classifier
# cnn_lowest = 0.50           # neurons with cnn probability lower than this value are rejected # adjusted


# opts_dict = {'fnames': fnames,
#             'fr': fr,
#             'decay_time': decay_time,
#             'strides': strides,
#             'overlaps': overlaps,
#             'max_shifts': max_shifts,
#             'max_deviation_rigid': max_deviation_rigid,
#             'pw_rigid': pw_rigid,
#             'p': p,
#             'nb': gnb,
#             'rf': rf,
#             'K': K, 
#             'stride': stride_cnmf,
#             'method_init': method_init,
#             'rolling_sum': True,
#             'only_init': True,
#             'ssub': ssub,
#             'tsub': tsub,
#             'merge_thr': merge_thr, 
#             'min_SNR': min_SNR,
#             'rval_thr': rval_thr,
#             'use_cnn': True,
#             'min_cnn_thr': cnn_thr,
#             'cnn_lowest': cnn_lowest}

# opts = params.CNMFParams(params_dict=opts_dict)


# In[5]:


#%% start a cluster for parallel processing (if a cluster already exists it will be closed and a new session will be opened)
if 'dview' in locals():
    cm.stop_server(dview=dview)
c, dview, n_processes = cm.cluster.setup_cluster(
    backend='local', n_processes=None, single_thread=False)


# ## Motion Correction & Memory mapping 
# 
# The cell below memory maps the file in order `'C'` and then loads the new memory mapped file. The saved files from motion correction are memory mapped files stored in `'F'` order. Their paths are stored in `mc.mmap_file`.

# In[6]:


mc = MotionCorrect(fnames, dview=dview, **opts.get_group('motion'))


# In[7]:


get_ipython().run_cell_magic('capture', '', '#%% Run piecewise-rigid motion correction using NoRMCorre\nmc.motion_correct(save_movie=True)\n')


# In[10]:


1281/60


# In[8]:


m_els = cm.load(mc.fname_tot_els)
border_to_0 = 0 if mc.border_nan is 'copy' else mc.border_to_0 
    # maximum shift to be used for trimming against NaNs

# 21h 34m - 23h 25m for 70k frame tif
# executed in 26h 13m (100k frame tiff)


# In[9]:


# fname_new = 'C:/Users/ll357/Documents/CaImAn/demos/temp_data/i1350_220225\\i1350_220225_002_multipage__d1_264_d2_796_d3_1_order_C_frames_80000_.mmap'


# In[10]:


# mc.mmap_file, fnames


# In[11]:


fname_new = cm.save_memmap(mc.mmap_file, base_name=fnames.split('\\')[-1][:-4] + '_', order='C',
                           border_to_0=border_to_0, dview=dview) # exclude borders
Yr, dims, T = cm.load_memmap(fname_new)
images = np.reshape(Yr.T, [T] + list(dims), order='F') 
    #load frames in python format (T x X x Y)
    
#%% restart cluster to clean up memory
cm.stop_server(dview=dview)
c, dview, n_processes = cm.cluster.setup_cluster(
    backend='local', n_processes=None, single_thread=False)


# ## Run CNMF on patches in parallel
# 
# - The FOV is split is different overlapping patches that are subsequently processed in parallel by the CNMF algorithm (nonnegative matrix factorization).
# - The results from all the patches are merged with special attention to idendtified components on the border.
# - The results are then refined by additional CNMF iterations.

# In[12]:


get_ipython().run_cell_magic('capture', '', 'cnm = cnmf.CNMF(n_processes, params=opts, dview=dview)\ncnm = cnm.fit(images)\n\n#%% plot contours of found components\nCn = cm.local_correlations(images.transpose(1,2,0))\nCn[np.isnan(Cn)] = 0\ncnm.estimates.plot_contours_nb(img=Cn)\n\n# 1h for 70k, 65-90 min for 80k frame tif\n')


# In[13]:


import pickle
def corr_img_save(file_now):
    Cn = cm.local_correlations(images.transpose(1,2,0))
    f = open(file_now[:-4] + '_Cn.pckl', 'wb')
    pickle.dump(Cn, f)
    f.close()
corr_img_save(fnames)
# todo: save correlation img in result .hdf5


# ## Component Evaluation
# 
# The processing in patches creates several spurious components. These are filtered out by evaluating each component using three different criteria:
# 
# - the shape of each component must be correlated with the data at the corresponding location within the FOV
# - a minimum peak SNR is required over the length of a transient
# - each shape passes a CNN based classifier

# In[14]:


cnm.estimates.Cn = Cn
cnm2 = cnm.refit(images, dview=dview) # cnm2 = cnm
cnm2.estimates.evaluate_components(images, cnm2.params, dview=dview)
print(len(cnm2.estimates.idx_components))
print(len(cnm2.estimates.idx_components_bad))
cnm2.estimates.plot_contours_nb(img=Cn, idx=cnm2.estimates.idx_components)


# In[15]:


plt.subplot(2,1,1)
plot_contours(cnm.estimates.A[:, cnm2.estimates.idx_components], Cn, 
              thr=0.9, display_numbers=False, 
              contour_args={'linewidth': 0.5, 'alpha': 0.5});
plt.subplot(2,1,2)
plot_contours(cnm.estimates.A[:, cnm2.estimates.idx_components_bad], Cn, 
              thr=0.9, display_numbers=False, 
              contour_args={'linewidth': 0.5, 'alpha': 0.5});
plt.savefig('Cn_' + fnames.split('\\')[-1][:-4] + '.pdf', bbox_inches='tight')


# In[16]:


#%% Extract DF/F values
cnm2.estimates.detrend_df_f(quantileMin=8, frames_window=250) # this step generates cnm2.F_dff
cnm2.estimates.F_dff.shape


# ## Closing, saving, and creating denoised version

# In[17]:


cnm2.save(fnames[:-4] + '.hdf5')

cm.stop_server(dview=dview)
log_files = glob.glob('*_LOG_*')
print(log_files)


# In[18]:


movie = cm.base.movies.movie(images,start_time=0,fr=30)
max_proj = np.max(movie, axis=0)
mean_img = np.mean(movie, axis=0)

plt.subplot(2, 1, 1)
plt.imshow(max_proj) # after registering
plt.subplot(2, 1, 2)
plt.imshow(mean_img)


# In[ ]:




