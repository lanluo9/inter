#!/usr/bin/env python
# coding: utf-8

# In[1]:

import os

img_dir = 'C:/Users/lan/Documents/repos/inter/code/py_playground/test_img' # absolute path of image folder
_, _, files = next(os.walk(img_dir)) # get files in directory
file_count = len(files) # file count is not recursive (does not count subdirectory)
print(file_count) # comment this to suppress output

# In[ ]:




