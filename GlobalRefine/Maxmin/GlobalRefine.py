# -*- coding: utf-8 -*-
"""
Created on Sun Apr  4 17:34:47 2021

@author: z5wu
"""

import os

import numpy as np
import matplotlib.pyplot as plt

import chaospy

import pickle

import scipy as sp

import random

#%% Determine the scriptpath
scriptPath = os.path.dirname(os.path.abspath(__file__)) + "\\"
os.chdir(scriptPath)

#%% Initial selected space
X_selected0 = pickle.load(open('TrainingSpace.sav', 'rb'))
# FillingPoints = pickle.load(open('FillingPoints.sav', 'rb'))

#%% Normalize the initial selected space
N_input0=len(X_selected0[0]) # Number of input variables
N_points0=len(X_selected0) # Number of training points

Selected0_normalize=np.zeros(shape=(N_points0,N_input0))
for i in range(N_input0):
    Input0_temp=X_selected0[:,i]
    Input0_max=np.amax(Input0_temp)
    Input0_min=np.amin(Input0_temp)
    Input0_normalize=(Input0_temp-Input0_min)/(Input0_max-Input0_min)
    Selected0_normalize[:,i]=Input0_normalize
    
#%% Downsample the current selected space
# N2use = 490 # Only use 500 out of all the 600 points
# random.seed(612)
# Index2use=random.sample(range(1, N_points0), 10) # Initial training points;

# for i in range(N2use):
#     X_initial=Selected0_normalize[Index2use,:] # Initial selected training points
#     dis_mat=sp.spatial.distance_matrix(X_initial,Selected0_normalize) # Distance matrix
#     dis_min=np.amin(dis_mat, 0)
#     index_new=np.argmax(dis_min)
#     Index2use.append(index_new)

# Selected0_normalize = Selected0_normalize[Index2use,:]

# #%% Save the selected training and testing space
# filename = 'Index2use.sav'
# pickle.dump(Index2use, open(filename, 'wb'))

#%% Create MC samples: Xtrain
n_samples = 1e5
distribution = chaospy.J(chaospy.Uniform(0,1),chaospy.Uniform(0,1),chaospy.Uniform(0,1),chaospy.Uniform(0,1),chaospy.Uniform(0,1))
Xtrain = distribution.sample(n_samples,rule="sobol")
Xtrain = np.transpose(Xtrain)

# Save the total MC points
filename = 'MC_space.sav'
pickle.dump(Xtrain, open(filename, 'wb'))


#%% Space Filling Parameters
N_input=len(Xtrain[0]) # Number of input variables
N_points=len(Xtrain) # Number of training points

N_select = 300 # The numbers of filling points
    
Index_selected=[] # Initial selected index;

#%%
for i in range(N_select):
    X_selected=np.append(Selected0_normalize,Xtrain[Index_selected,:],axis=0) # Selected training points
    dis_mat=sp.spatial.distance_matrix(X_selected,Xtrain) # Distance matrix
    dis_min=np.amin(dis_mat, 0)
    index_new=np.argmax(dis_min)
    Index_selected.append(index_new)
    
#%% Denormalization
FillingPoints = np.zeros(shape=(len(Index_selected),N_input0))
for i in range(N_input0):
    Input_temp=Xtrain[Index_selected,i]
    Input0_temp=X_selected0[:,i]
    Input0_max=np.amax(Input0_temp)
    Input0_min=np.amin(Input0_temp)
    Input_denormalize=Input_temp*(Input0_max-Input0_min)+Input0_min
    FillingPoints[:,i]=Input_denormalize
    

#%% Save the filling points
filename = 'FillingPoints300.sav'
pickle.dump(FillingPoints, open(filename, 'wb'))

#%% Save the denormed MC space
Xtrain_denorm = np.zeros(shape=(N_points,N_input0))
for i in range(N_input0):
    Input_temp=Xtrain[:,i]
    Input0_temp=X_selected0[:,i]
    Input0_max=np.amax(Input0_temp)
    Input0_min=np.amin(Input0_temp)
    Input_denormalize=Input_temp*(Input0_max-Input0_min)+Input0_min
    Xtrain_denorm[:,i]=Input_denormalize

filename = 'MC_denorm.sav'
pickle.dump(Xtrain_denorm, open(filename, 'wb'))
    
    
#%% Visualization
# Selected0_normalize = X_selected0[Index2use,:]
# X1 = Selected0_normalize[:,0]
# X2 = Selected0_normalize[:,1]
# X3 = Selected0_normalize[:,2]
# X4 = Selected0_normalize[:,3]
# X5 = Selected0_normalize[:,4]
# X6 = Selected0_normalize[:,5]
# X1new = FillingPoints[:,0]
# X2new = FillingPoints[:,1]
# X3new = FillingPoints[:,2]
# X4new = FillingPoints[:,3]
# X5new = FillingPoints[:,4]
# X6new = FillingPoints[:,5]
# #%% 1&2
# fig, ax = plt.subplots()
# ax.scatter(X1, X2,label='Training Points',s=5)
# ax.scatter(X1new, X2new, label='Added Points', color='red',s=20,marker='x')
# ax.legend()
# plt.title('Dimension 1&2')
# plt.xlabel('Dimension 1')
# plt.ylabel('Dimension 2')
# plt.show()

#%% 1&3
# fig, ax = plt.subplots()
# ax.scatter(X1, X3,label='Training Points',s=5)
# ax.scatter(X1new, X3new, label='Added Points', color='red',s=20,marker='x')
# ax.legend()
# plt.title('Dimension 1&3')
# plt.xlabel('Dimension 1')
# plt.ylabel('Dimension 3')
# plt.show()

# #%% 1&4
# fig, ax = plt.subplots()
# ax.scatter(X1, X4,label='Training Points',s=5)
# ax.scatter(X1new, X4new, label='Added Points', color='red',s=20,marker='x')
# ax.legend()
# plt.title('Dimension 1&4')
# plt.xlabel('Dimension 1')
# plt.ylabel('Dimension 4')
# plt.show()


# #%% 2&3
# fig, ax = plt.subplots()
# ax.scatter(X2, X3,label='Training Points',s=5)
# ax.scatter(X2new, X3new, label='Added Points', color='red',s=20,marker='x')
# ax.legend()
# plt.title('Dimension 2&3')
# plt.xlabel('Dimension 2')
# plt.ylabel('Dimension 3')
# plt.show()



