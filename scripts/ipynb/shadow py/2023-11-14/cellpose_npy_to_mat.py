#!/usr/bin/env python
# coding: utf-8

# In[28]:


import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio

dir_sess = r'Z:\All_Staff\home\lan\Analysis\2P\210113_i1329\210113_i1329_runs-006'.replace('\\', '/')
npy_list = [fn for fn in os.listdir(dir_sess) if fn.endswith('.npy')]
npy_list = [fn for fn in npy_list if 'v1' not in fn] # drop item if contains 'v1', which is cellpose version 1, auto seg. bad for LI data
assert len(npy_list) == 1

cellpose_npy = npy_list[0]
npy = np.load(os.path.join(dir_sess, cellpose_npy), allow_pickle=True)
# npy.item().keys()
# npy.item()['masks'].shape, npy.item()['ismanual']

plt.imshow(npy.item()['masks'])

# save as mat
masks = npy.item()['masks']
assert masks.shape == (264, 796)
sio.savemat(os.path.join(dir_sess, 'cellpose_mask.mat'), mdict={'cellpose_mask': masks}) # save cellpose mask to mat

