# -*- coding: utf-8 -*-
"""
Created on Wed May 19 14:33:56 2021

@author: MDT
"""
import os
import numpy as np
import matplotlib.pyplot as plt

import pickle

#%% Determine the scriptpath
scriptPath = os.path.dirname(os.path.abspath(__file__)) + "\\"
os.chdir(scriptPath)

workPath = scriptPath

        
#%% Extracting each MSE
MSE = []

N_max = 5

for i in range(0,N_max):
    MSE_folder = workPath+"Model{}\\".format(i)
    os.chdir(MSE_folder)
    MSE_temp = pickle.load(open('Var_max.sav', 'rb'))
    MSE.append(MSE_temp)
    
#%%
aaa = MSE
plt.figure()
plt.plot(np.array(range(0,N_max)),MSE)
plt.xlabel('Number of added points')
plt.ylabel('$MSE$')