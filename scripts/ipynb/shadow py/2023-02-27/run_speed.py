#!/usr/bin/env python
# coding: utf-8

# In[9]:


import numpy as np
import scipy.io as sio


# In[5]:


mouse = 'i1375'
date = '220915'
run = '003'
time = '1421'

data_file = r'\\duhs-user-nc1.dhe.duke.edu\dusom_glickfeldlab\All_staff\Behavior\Data\data-' + mouse + '-' + date + '-' + time + '.mat'
data_file


# In[11]:


data = sio.loadmat(data_file)
# data['input']

# # Load the data into Python     
# D= sio.loadmat('data.mat')

# build a list of keys and values for each entry in the structure
vals = data['input'][0,0] #<-- set the array you want to access. 
keys = data['input'][0,0].dtype.descr

# Assemble the keys and values into variables with the same name as that used in MATLAB
for i in range(len(keys)):
    key = keys[i][0]
    # val = np.squeeze(vals[key][0][0])  # squeeze is used to covert matlat (1,n) arrays into numpy (1,) arrays. 
    exec(key + '=val')


# In[13]:


def _check_keys( dict):
    """
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    """
    for key in dict:
        if isinstance(dict[key], sio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict


def _todict(matobj):
    """
    A recursive function which constructs from matobjects nested dictionaries
    """
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, sio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict


def loadmat(filename):
    """
    this function should be called instead of direct scipy.io .loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    """
    data = sio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)


# In[15]:


data = loadmat(data_file)
data['input']

