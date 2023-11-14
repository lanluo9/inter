#!/usr/bin/env python
# coding: utf-8

# In[1]:


from PIL import Image, ImageDraw, ImageOps
import cv2 

import os
from os import listdir
from os.path import isfile, join
import glob
from tqdm.notebook import tqdm

import numpy as np
import math
import matplotlib.pyplot as plt
import seaborn as sns


# # unify lum dist
# to unify mean overall for image set, 
# check population mean -> calc indiv img mean -> shift img mean towards pop mean

# In[56]:


os.chdir(r'Z:\All_Staff\home\lan\Mwork\mix50 - bunnytop high lum contrast mix grating and noise\Image'.replace('\\','/'))
img = Image.open('1.png')
im_arr = np.array(img)
plt.imshow(im_arr); plt.colorbar()


# In[57]:


indiv_img_mean = np.mean(im_arr)
indiv_img_mean


# In[58]:


# imagenet stats from torchvision.transforms.Normalize()
channel_mean = np.array([0.485, 0.456, 0.406]) * 255 # pytorch tensor range 0-1 convert to pixel value 0-255
pop_img_mean = np.mean(channel_mean)
pop_img_mean


# In[59]:


im_arr[:,:,0:3] = im_arr[:,:,0:3] + np.uint8(pop_img_mean - indiv_img_mean)
im_arr.shape
plt.imshow(im_arr); plt.colorbar()

# # TODO: after shifting mean, prevent out of range values get clipped


# # quantify coverage of feature space (esp spatial freq) in img set
# sum of all img SF power spectrum -> test uniform distribution -> plot in order  
# single img SF power spectrum -> test unimodal distribution -> test peak -> plot in order

# In[ ]:




