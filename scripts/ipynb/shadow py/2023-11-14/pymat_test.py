#!/usr/bin/env python
# coding: utf-8

# # test pymat bridge

# In[1]:


from pymatbridge import Matlab
mlab = Matlab()


# In[2]:


mlab.start()


# In[3]:


results = mlab.run_code('a=1')
# print(results)
mlab.get_variable('a')


# In[4]:


dir_func = r'C:\Users\GlickfeldLab\Documents\test\inter\scripts\misc\personalized_matlab\draft.m'.replace('\\', '/')
res = mlab.run_func(dir_func, {'arg1': 3, 'arg2': 5})
print(res['result'])
res


# In[5]:


mlab.stop()


# # cell magic bug

# In[7]:


get_ipython().run_line_magic('load_ext', 'pymatbridge')


# In[8]:


get_ipython().run_cell_magic('matlab', '', "\ndisp('hello');\n\na = linspace(0.01,6*pi,100);\nplot(sin(a))\ngrid on\nhold on\nplot(cos(a),'r')\n")

