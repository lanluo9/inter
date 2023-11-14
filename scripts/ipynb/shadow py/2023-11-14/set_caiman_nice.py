#!/usr/bin/env python
# coding: utf-8

# In[2]:


import psutil
import time
import os
from IPython.display import clear_output

clear_output(wait=False)

current_pid = os.getpid()
print(f'setting high prio for set_caiman_nice self pid = {current_pid}')
psutil.Process(os.getpid()).nice(psutil.ABOVE_NORMAL_PRIORITY_CLASS)

local_py_pids = []
other_user_pids = []
while True:
# while psutil.cpu_percent() > 0:
    print(f'cpu util = {psutil.cpu_percent()}, start setting low prio for current user\'s caiman process.')

    for proc in psutil.process_iter(): # reading psutil.process_iter into a df.DataFrame is even slower
        try:
            if ("python" in proc.name()) and (os.getlogin() in proc.username() and (proc.pid not in other_user_pids)):
                local_py_pids.append(proc.pid)
        except psutil.AccessDenied:
            other_user_pids.append(proc.pid)
            pass

    local_py_pids.remove(os.getpid()) # get pid for all current user's python processes except for this script itself
    for pid in local_py_pids:
        if psutil.Process(pid).nice() != psutil.BELOW_NORMAL_PRIORITY_CLASS: # check prio is faster than setting prio
            psutil.Process(pid).nice(psutil.BELOW_NORMAL_PRIORITY_CLASS)
            print("lower prio of process {} owned by {}".format(pid, psutil.Process(pid).username()))


    time.sleep(60*60) # check for new python process every n sec
    clear_output(wait=True) # clear long cell output after pause

