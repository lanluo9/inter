#!/usr/bin/env python
# coding: utf-8

# # objective
# check cellpose use case for multiple image / multisess alignment (kicking out repeated cells)

# In[1]:


import matplotlib.pyplot as plt
import numpy as np
import os
# from skimage.segmentation import find_boundaries


# # try w matlab register ref
# ## test merge ROI

# In[64]:


import scipy.io as sio
import os

os.chdir(r'Z:\All_Staff\home\lan\Analysis\2P\220711_i1372\220711_i1372_runs-002'.replace('\\','/'))
data = sio.loadmat('cellpose_mask.mat')
cellpose_mask_0 = data['cellpose_mask']

os.chdir(r'Z:\All_Staff\home\lan\Analysis\2P\220714_i1372\220714_i1372_runs-003'.replace('\\','/'))
data = sio.loadmat('cellpose_mask.mat')
cellpose_mask_1 = data['cellpose_mask']
cellpose_mask_1.shape


# In[66]:


# modify mask 0 and 1 to make mask sum have fewer cells

# cellpose_mask_0[40:, 50:] = 0
cellpose_mask_0[cellpose_mask_0 > 10] = 0
cellpose_mask_1[:-40, :-50] = 0
plt.imshow(np.array(cellpose_mask_0, dtype=bool));
plt.imshow(np.array(cellpose_mask_1, dtype=bool));


# In[68]:


plt.figure(figsize=(20,8))
mask_sum = cellpose_mask_0 + cellpose_mask_1
plt.imshow(np.array(mask_sum, dtype=bool))
plt.colorbar()


# In[71]:


import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

from skimage import data
from skimage.filters import threshold_otsu
from skimage.segmentation import clear_border
from skimage.measure import label, regionprops
from skimage.morphology import closing, square
from skimage.color import label2rgb

data = np.array(mask_sum, dtype=bool) #*255.
# data = np.uint8(data)
image = data # data.coins()[50:-50, 50:-50]

# apply threshold
thresh = threshold_otsu(image)
bw = closing(image > thresh, square(3))

# remove artifacts connected to image border
cleared = clear_border(bw)

# label image regions
label_image = label(cleared)
# to make the background transparent, pass the value of `bg_label`,
# and leave `bg_color` as `None` and `kind` as `overlay`
image_label_overlay = label2rgb(label_image, image=image, bg_label=0)

fig, ax = plt.subplots(figsize=(10, 6))
ax.imshow(image_label_overlay)
# ax.colorbar()

# for region in regionprops(label_image):
#     # take regions with large enough areas
#     if region.area >= 100:
#         # draw rectangle around segmented coins
#         minr, minc, maxr, maxc = region.bbox
#         rect = mpatches.Rectangle((minc, minr), maxc - minc, maxr - minr,
#                                   fill=False, edgecolor='red', linewidth=2)
#         ax.add_patch(rect)

# ax.set_axis_off()
plt.tight_layout()
plt.show()


# In[75]:


plt.figure(figsize=(20,8))
plt.imshow(label_image, cmap='plasma')
plt.colorbar();

np.unique(label_image)


# # try w FAIMCalcium

# In[1]:


os.chdir(r'C:\Users\lan\Documents\repos\inter\scripts\misc\borrowed\FAIMCalcium-master'.replace('\\','/'))
# conda install -c menpo opencv
from AffineCa2p.FAIM import FAIMCaSig

import cv2
cv2.xfeatures2d


# ## tutorial run

# In[4]:


dir_A5 = r'C:\Users\lan\Documents\repos\inter\scripts\misc\borrowed\FAIMCalcium-master\examples\A5'.replace('\\','/')
Tmatrices, regImages, regROIs = FAIMCaSig.AlignIm(dir_A5)

# errored out on 11 min
# ASIFT is running
# registering reg_day08A5.png
# FAIM\find_obj.py in init_feature(name)
# module 'cv2' has no attribute 'xfeatures2d'


# ## result

# In[11]:


len(Tmatrices), Tmatrices[0].shape, len(regImages), regImages[0].shape, len(regROIs), regROIs[0].shape


# In[15]:


plt.subplot(1,2,1)
plt.imshow(regROIs[0])
plt.subplot(1,2,2)
plt.imshow(regROIs[1])


# In[28]:


output_Im = np.zeros([np.size(regROIs[0],0), np.size(regROIs[0],1), 3], np.uint8)

outlines1 = regROIs[0]
# outlines1 = np.zeros(regROIs[0].shape, bool)
# outlines1[find_boundaries(regROIs[0], mode='inner')] = 1
outX1, outY1 = np.nonzero(outlines1)
output_Im[outX1, outY1] = np.array([255, 0, 0])

outlines2 = regROIs[1]
# outlines2 = np.zeros(regROIs[1].shape, bool)
# outlines2[find_boundaries(regROIs[1], mode='inner')] = 1
outX2, outY2 = np.nonzero(outlines2)
output_Im[outX2, outY2] = np.array([255, 255, 22])

plt.imshow(output_Im)


# In[16]:


import matplotlib.pyplot as plt
plt.subplot(1,2,1)
plt.imshow(regImages[0])
plt.subplot(1,2,2)
plt.imshow(regImages[1])


# In[ ]:




