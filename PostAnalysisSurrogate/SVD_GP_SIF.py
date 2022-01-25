# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 10:19:29 2021

@author: MDT
"""
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel as C_ker

import chaospy

from sklearn.metrics import mean_squared_error

import pickle

import random

random.seed(612)
np.random.seed(612)

#%% Determine the scriptpath
scriptPath = os.path.dirname(os.path.abspath(__file__)) + "\\"
os.chdir(scriptPath)

displPath = scriptPath + "LocalAppliedDispl\\"

SIFPath = scriptPath + "SIF\\"

#%% Read all existing indexes from SIF folder

N_total = 200

indexList = []


OutputSpace = []

for index in range(0,N_total):
    fileName = SIFPath + "S" + str(index) + ".csv"
    if os.path.exists(fileName):
        indexList += [index]
        SIF_full = np.array(pd.read_csv(fileName))
        OutputSpace = np.append(OutputSpace,np.mean(SIF_full[:,0]))

#%% Read the uniform pressure and the crack length
n_samples = 200

distribution = chaospy.Iid(chaospy.Uniform(0, 1), 4)
samples = distribution.sample(n_samples,rule="sobol")
Input4 = 0.5+samples[3,:]*3.5
        
#%% Extracting the last iteration of each LocalAppliedDispl as input
InputSpace = np.zeros((120*6+1,1))
max_interation = 15

for i in indexList:
    j = 0
    InputName = displPath+"S{}_Displ{}.csv".format(i,(j))
    while os.path.exists(InputName):
        j = j+1
        InputName = displPath+"S{}_Displ{}.csv".format(i,(j))
    else:
        j = j-1
        InputName = displPath+"S{}_Displ{}.csv".format(i,(j))
        U_index = pd.read_csv(InputName,usecols = ['0','1','2','3','4','5'])
        U_array = np.reshape(np.array(U_index),(120*6,1))
        U_array = np.append(U_array,Input4[i,])
        InputSpace = np.append(InputSpace,np.reshape(U_array,(120*6+1,1)),axis=1)
          
InputSpace = np.transpose(InputSpace[:,1:])

#%% Define train, test and vali
N_train = 185
N_test = 10
N_vali = len(OutputSpace)-N_train-N_test

index_total=random.sample(range(0, len(OutputSpace)), N_train+N_test)
Input_total = InputSpace[index_total,:]
Output_total = OutputSpace[index_total]

index_vali = np.delete(range(0, len(OutputSpace)), index_total)
Input_vali = InputSpace[index_vali,:]
Output_vali = OutputSpace[index_vali]


#%% SVD for both input datasets
U1,S1,V1 =  np.linalg.svd(Input_total[:,0:120*6],full_matrices=True)

#%% SVD for input dataset
vars1 = np.diag(S1)

totalVar1 = np.cumsum(S1)/np.sum(S1)
Rank1 = 4
#%%
Reduced_Input_train_0 = U1[0:N_train, :Rank1]
Reduced_Input_train = np.append(Reduced_Input_train_0,np.reshape(Input_total[0:N_train,120*6],(len(Reduced_Input_train_0[:,1]),1)),axis=1)

Reduced_Input_test_0 = U1[N_train:, :Rank1]
Reduced_Input_test = np.append(Reduced_Input_test_0,np.reshape(Input_total[N_train:,120*6],(len(Reduced_Input_test_0[:,1]),1)),axis=1)

encoder1 = np.linalg.pinv(np.diag(S1[:Rank1]) @ V1[:Rank1, :])
filename = 'InputEncoder.sav'
pickle.dump(encoder1, open(filename, 'wb'))
decoder1 = np.diag(S1[:Rank1]) @ V1[:Rank1, :]
filename = 'InputDecoder.sav'
pickle.dump(decoder1, open(filename, 'wb'))


#%% Train the surrogate model
nInputs=len(Reduced_Input_train[0]) # Number of input variables
samples_i = np.array(range(N_train+N_test+N_vali))

Output_train = Output_total[:N_train]
Output_test = Output_total[N_train:]

#%% Surrogate 1 for Output1
L0=10
L_bounds=(1e-5,1e5)

trainMSE_1=[]  
testMSE_1=[] 

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

gp_GPML1 = GaussianProcessRegressor(kernel=kernel_prod, alpha=1e-7,n_restarts_optimizer=100,normalize_y=False,random_state=612)
gp_GPML1.fit(Reduced_Input_train, Output_train)
     
y_output = gp_GPML1.predict(Reduced_Input_train, return_std=False)
trainMSE_1.append(np.sqrt(mean_squared_error(Output_train, y_output, multioutput='uniform_average')))

y_pred, sigma = gp_GPML1.predict(Reduced_Input_test, return_std=True)
testMSE_1.append(np.sqrt(mean_squared_error(Output_test, y_pred, multioutput='uniform_average')))

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(samples_i[:N_train],Output_train, label='True Output1')
plt.plot(samples_i[:N_train],y_output, label='Predicted Output1')
plt.xlabel('$Sample Points$')
plt.ylabel('$Output$')
# plt.xlim(120,180)
plt.legend(loc='upper right')
plt.subplot(2, 1, 2)
plt.plot(samples_i[N_train:-N_vali],Output_test, label='True Output1')
plt.plot(samples_i[N_train:-N_vali],y_pred, label='Predicted Output1')
plt.fill(np.concatenate([samples_i[N_train:-N_vali], samples_i[N_train:-N_vali][::-1]]),np.concatenate([y_pred - 3 * sigma, (y_pred + 3 * sigma)[::-1]]), alpha=.3, fc='b', ec='None')
plt.xlabel('$Sample Points$')
plt.ylabel('$Output$')
plt.legend(loc='upper right')
    
print("training error" + str(trainMSE_1))
print("testing error" + str(testMSE_1))


#%% Save all files in one variable
GP_list = [encoder1,decoder1,gp_GPML1]
file_name = "Post_SIF.pkl"

pickle.dump(GP_list, open(file_name, "wb"))


#%% Load the surrogate models
 
GP_Model = pickle.load(open("Post_SIF.pkl", "rb"))
encoder1 = GP_Model[0]
decoder1 = GP_Model[1]
Model1 = GP_Model[2]



#%% Validation

Reduced_Input_vali_0 = Input_vali[:,0:120*6] @ encoder1
Reduced_Input_vali = np.append(Reduced_Input_vali_0,np.reshape(Input_vali[:,120*6],(len(Reduced_Input_vali_0[:,1]),1)),axis=1)

Reduced_Output1 = Model1.predict(Reduced_Input_vali,return_std=False).reshape(N_vali,1)


#%% Test
testSpace = np.zeros((120*6+1,1))
for i in range(10):
    fileName = scriptPath+"S0_Displ{}.csv".format(i)
    U_index = pd.read_csv(fileName,usecols = ['0','1','2','3','4','5'])
    U_array = np.reshape(np.array(U_index),(120*6,1))
    U_array = np.append(U_array,1)
    testSpace = np.append(testSpace,np.reshape(U_array,(120*6+1,1)),axis=1)


testSpace = np.transpose(testSpace[:,1:])
#%%
Reduced_Input_vali_0 = testSpace[:,0:120*6] @ encoder1
Reduced_Input_vali = np.append(Reduced_Input_vali_0,np.reshape(testSpace[:,120*6],(len(Reduced_Input_vali_0[:,1]),1)),axis=1)

Reduced_Output1 = Model1.predict(Reduced_Input_vali,return_std=False).reshape(10,1)



