# -*- coding: utf-8 -*-
"""
Created on Mon Apr 12 16:20:32 2021

@author: MDT
"""

import numpy as np
import pickle
import os
import sys

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel as C_ker


#%% The passing variables
index = int(sys.argv[-2])
# index = 0

scriptPath = sys.argv[-1]
# scriptPath = "C:\ZIHAN\GreenupCrackEstimation\A1201\LocalRefine\MSE\\"
os.chdir(scriptPath)

#%% Set the scriptpath and workpath

workPath= scriptPath + "\\Files\\"

modelPath = workPath + "\\Model" + str(index)

#%% Load current sample points
os.chdir(modelPath)
Reduced_Input_train = pickle.load(open('Input_train_Global200.sav', 'rb'))
Reduced_Output_train = pickle.load(open('Output_train_Global200.sav', 'rb'))

#%%
nInputs=len(Reduced_Input_train[0]) # Number of input variables

#%% Update the parameters
MSE1 = []
MSE2 = []
MSE3 = []
MSE4 = []
MSE5 = []


# Surrogate 1 for Output1
L0=10
L_bounds=(1e-5,1e5)

############### Constant Kernel##############
kernel_Con=C_ker(L0, (1e-5,1e5))
############### Matern Kernel##############
Matern_length0=[L0]
Matern_length_bounds=[L_bounds]
for i in range(nInputs-1):
    Matern_length0.append(L0)
    Matern_length_bounds.append(L_bounds)
kernel_Matern = Matern(Matern_length0, Matern_length_bounds, nu=2.5)
kernel_prod=kernel_Con*kernel_Matern

gp_GPML1 = GaussianProcessRegressor(kernel=kernel_prod, alpha=1e-8,n_restarts_optimizer=10,normalize_y=False,random_state=612)
gp_GPML1.fit(Reduced_Input_train, Reduced_Output_train[:,0])


# Surrogate 2 for Output2
L0=10
L_bounds=(1e-5,1e5)

############### Constant Kernel##############
kernel_Con=C_ker(L0, (1e-5,1e5))
############### Matern Kernel##############
Matern_length0=[L0]
Matern_length_bounds=[L_bounds]
for i in range(nInputs-1):
    Matern_length0.append(L0)
    Matern_length_bounds.append(L_bounds)
kernel_Matern = Matern(Matern_length0, Matern_length_bounds, nu=2.5)
kernel_prod=kernel_Con*kernel_Matern

gp_GPML2 = GaussianProcessRegressor(kernel=kernel_prod, alpha=1e-8,n_restarts_optimizer=10,normalize_y=False,random_state=123)
gp_GPML2.fit(Reduced_Input_train, Reduced_Output_train[:,1])


# Surrogate 3 for Output3
L0=10
L_bounds=(1e-5,1e5)

############### Constant Kernel##############
kernel_Con=C_ker(L0, (1e-5,1e5))
############### Matern Kernel##############
Matern_length0=[L0]
Matern_length_bounds=[L_bounds]
for i in range(nInputs-1):
    Matern_length0.append(L0)
    Matern_length_bounds.append(L_bounds)
kernel_Matern = Matern(Matern_length0, Matern_length_bounds, nu=1.5)
kernel_prod=kernel_Con*kernel_Matern

gp_GPML3 = GaussianProcessRegressor(kernel=kernel_prod, alpha=1e-8,n_restarts_optimizer=10,normalize_y=False,random_state=612)
gp_GPML3.fit(Reduced_Input_train, Reduced_Output_train[:,2])


# Surrogate 4 for Output4
L0=10
L_bounds=(1e-5,1e5)

############### Constant Kernel##############
kernel_Con=C_ker(L0, (1e-5,1e5))
############### Matern Kernel##############
Matern_length0=[L0]
Matern_length_bounds=[L_bounds]
for i in range(nInputs-1):
    Matern_length0.append(L0)
    Matern_length_bounds.append(L_bounds)
kernel_Matern = Matern(Matern_length0, Matern_length_bounds, nu=2.5)
kernel_prod=kernel_Con*kernel_Matern

gp_GPML4 = GaussianProcessRegressor(kernel=kernel_prod, alpha=1e-8,n_restarts_optimizer=10,normalize_y=False,random_state=612)
gp_GPML4.fit(Reduced_Input_train, Reduced_Output_train[:,3])


# Surrogate 5 for Output4
L0=10
L_bounds=(1e-5,1e5)

############### Constant Kernel##############
kernel_Con=C_ker(L0, (1e-5,1e5))
############### Matern Kernel##############
Matern_length0=[L0]
Matern_length_bounds=[L_bounds]
for i in range(nInputs-1):
    Matern_length0.append(L0)
    Matern_length_bounds.append(L_bounds)
kernel_Matern = Matern(Matern_length0, Matern_length_bounds, nu=2.5)
kernel_prod=kernel_Con*kernel_Matern

gp_GPML5 = GaussianProcessRegressor(kernel=kernel_prod, alpha=1e-8,n_restarts_optimizer=10,normalize_y=False,random_state=612)
gp_GPML5.fit(Reduced_Input_train, Reduced_Output_train[:,4])


#%% Obtain the vairance from the model
os.chdir(scriptPath)
Xtrain = pickle.load(open('MC_denorm.sav', 'rb'))

Output1,sigma1 = gp_GPML1.predict(Xtrain,return_std=True)
Output2,sigma2 = gp_GPML2.predict(Xtrain,return_std=True)
Output3,sigma3 = gp_GPML3.predict(Xtrain,return_std=True)
Output4,sigma4 = gp_GPML4.predict(Xtrain,return_std=True)
Output5,sigma5 = gp_GPML5.predict(Xtrain,return_std=True)
#%%
MeanVar = np.mean([sigma1,sigma2,sigma3,sigma4,sigma5], axis = 0)
os.chdir(modelPath)
filename = 'Var_max.sav'
Var_max = np.amax(MeanVar)
pickle.dump(Var_max, open(filename, 'wb'))
#%%
index_Updating = np.argmax(MeanVar)
AddedPoints = []        
AddedPoints.append(Xtrain[index_Updating,:])

#%%
nextModelPath = workPath + "\\Model" + str(index+1)
if not os.path.exists(nextModelPath):
    os.makedirs(nextModelPath) 
os.chdir(nextModelPath)
filename = 'AddedPoints.sav'
AddedPoints = np.array(AddedPoints)
pickle.dump(AddedPoints, open(filename, 'wb'))