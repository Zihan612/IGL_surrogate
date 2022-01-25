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

reactionPath = scriptPath + "LocalReaction\\"

displPath = scriptPath + "LocalAppliedDispl\\"

SIFPath = scriptPath + "SIF\\"

#%% Read all existing indexes from SIF folder

N_total = 200

indexList = []

for index in range(0,N_total):
    fileName = SIFPath + "S" + str(index) + ".csv"
    if os.path.exists(fileName):
        indexList += [index]

#%% Read the uniform pressure and the crack length
n_samples = 200

distribution = chaospy.Iid(chaospy.Uniform(0, 1), 4)
samples = distribution.sample(n_samples,rule="sobol")

# Hydraulic pressure h_up
Input1 = 432+samples[0,:]*288 # range:36*12-60*12

# Hydraulic pressure h_down
Input2 = 120+samples[1,:]*240 # range:10*12-30*12

# Gap length
Input3 = samples[2,:]*150

# Crack length
Input4 = 0.5+samples[3,:]*3.5
        
#%% Extracting all the LocalReaction and LocalAppliedDispl as output and input, respectively

nDOF = 120*6 # Total DOFs

U_total = np.zeros((nDOF+1,1))
P_total = np.zeros((nDOF,1))

max_interation = 10

for i in indexList:
    for j in range(0,max_interation):
        OutputName = reactionPath+"S{}_Reaction{}.csv".format(i,j)
        InputName = displPath+"S{}_Displ{}.csv".format(i,(j))
        if os.path.exists(OutputName) and os.path.exists(OutputName):
            
            P_index = pd.read_csv(OutputName,usecols = ['0','1','2','3','4','5'])
            P_array = np.reshape(np.array(P_index),(nDOF,1))
            P_total = np.append(P_total,P_array,axis=1)
            
            U_index = pd.read_csv(InputName,usecols = ['0','1','2','3','4','5'])
            U_array = np.reshape(np.array(U_index),(nDOF,1))
            U_array = np.append(U_array,Input4[i,])
            U_total = np.append(U_total,np.reshape(U_array,(nDOF+1,1)),axis=1)
            
P_total = np.transpose(P_total[:,1:])
U_total = np.transpose(U_total[:,1:])

#%% Seperating the dataset into training, testing, and validating
N_all = len(P_total)

N_train = 700
N_test = 30
N_vali = N_all-N_train-N_test

index_train=random.sample(range(0, N_all), N_train+N_test)
input_train = U_total[index_train,:]
output_train = P_total[index_train,:]

index_vali = np.delete(range(0, N_all), index_train)
input_vali = U_total[index_vali,:]
output_vali = P_total[index_vali,:]


#%% SVD for both input and output datasets
U1,S1,V1 =  np.linalg.svd(input_train[:,:-1],full_matrices=True)
U2,S2,V2 =  np.linalg.svd(output_train,full_matrices=True)

#%% SVD for input dataset
vars1 = np.diag(S1)

totalVar1 = np.cumsum(S1)/np.sum(S1)
Rank1 = 4
#%%
Reduced_Input_train_0 = U1[0:N_train, :Rank1]
Reduced_Input_train = np.append(Reduced_Input_train_0,np.reshape(input_train[0:N_train,nDOF],(len(Reduced_Input_train_0[:,1]),1)),axis=1)

Reduced_Input_test_0 = U1[N_train:, :Rank1]
Reduced_Input_test = np.append(Reduced_Input_test_0,np.reshape(input_train[N_train:,nDOF],(len(Reduced_Input_test_0[:,1]),1)),axis=1)

encoder1 = np.linalg.pinv(np.diag(S1[:Rank1]) @ V1[:Rank1, :])
# filename = 'InputEncoder.sav'
# pickle.dump(encoder1, open(filename, 'wb'))
decoder1 = np.diag(S1[:Rank1]) @ V1[:Rank1, :]
# filename = 'InputDecoder.sav'
# pickle.dump(decoder1, open(filename, 'wb'))

#%% SVD for output dataset
vars2 = np.diag(S2)

totalVar2 = np.cumsum(S2)/np.sum(S2)

Rank2 = 5
#%%
Reduced_Output_train = U2[0:N_train, :Rank2]
Reduced_Output_test = U2[N_train:, :Rank2]

decoder2 = np.diag(S2[:Rank2]) @ V2[:Rank2, :]
# filename = 'OutputDecoder.sav'
# pickle.dump(decoder2, open(filename, 'wb'))
encoder2 = np.linalg.pinv(decoder2)
# filename = 'OutputEncoder.sav'
# pickle.dump(encoder2, open(filename, 'wb'))

#%% Add New Training Samples
#%% Read all existing indexes from DataSaved
NewOutputPath = scriptPath + "DataSaved\\"

N_total_New = 300

indexListNew = []

for index in range(0,N_total_New):
    fileName = NewOutputPath + "S" + str(index) + "_Reaction.csv"
    if os.path.exists(fileName):
        indexListNew += [index]

#%% Read New Input Space
InputSpace = pickle.load(open('FillingPoints300.sav', 'rb'))

#%%
InputUsed = InputSpace[indexListNew,:]

# Index2use = pickle.load(open('Index2use.sav', 'rb'))

# Reduced_Input_train = Reduced_Input_train[Index2use,:]

Reduced_Input_train = np.append(Reduced_Input_train,InputUsed,axis=0)
                                
#%% Read New Output Space
P_total = np.zeros((nDOF,1))

for i in indexListNew:
    OutputName = NewOutputPath+"S{}_Reaction.csv".format(i)        
    P_index = pd.read_csv(OutputName,usecols = ['0','1','2','3','4','5'])
    P_array = np.reshape(np.array(P_index),(nDOF,1))
    P_total = np.append(P_total,P_array,axis=1)
                   
P_total = np.transpose(P_total[:,1:])

Reduced_Output_New = P_total @ encoder2

# Reduced_Output_train = Reduced_Output_train[Index2use,:]

Reduced_Output_train = np.append(Reduced_Output_train,Reduced_Output_New,axis=0)

#%% Determine the number of global refine point to use
filename = 'Input_trainNew200.sav'
Reduced_Input_train=Reduced_Input_train[:-100,:]
pickle.dump(Reduced_Input_train, open(filename, 'wb'))
filename = 'Output_trainNew200.sav'
Reduced_Output_train=Reduced_Output_train[:-100,:]
pickle.dump(Reduced_Output_train, open(filename, 'wb'))


#%% Train the surrogate model
nInputs=len(Reduced_Input_train[0]) # Number of input variables

#%% Update the parameters
samples_i = np.array(range(len(Reduced_Output_train[:,0])+N_test+N_vali))
N_train = len(Reduced_Output_train[:,0])

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

gp_GPML1 = GaussianProcessRegressor(kernel=kernel_prod, alpha=1e-8,n_restarts_optimizer=0,normalize_y=False,random_state=612)
gp_GPML1.fit(Reduced_Input_train, Reduced_Output_train[:,0])
     
y_output = gp_GPML1.predict(Reduced_Input_train, return_std=False)
trainMSE_1.append(np.sqrt(mean_squared_error(Reduced_Output_train[:,0], y_output, multioutput='uniform_average')))

y_pred, sigma = gp_GPML1.predict(Reduced_Input_test, return_std=True)
testMSE_1.append(np.sqrt(mean_squared_error(Reduced_Output_test[:,0], y_pred, multioutput='uniform_average')))

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(samples_i[:N_train],Reduced_Output_train[:,0], label='True Output1')
plt.plot(samples_i[:N_train],y_output, label='Predicted Output1')
plt.xlabel('$Sample Points$')
plt.ylabel('$Output$')
# plt.xlim(1000,1050)
plt.legend(loc='upper right')
plt.subplot(2, 1, 2)
plt.plot(samples_i[N_train:-N_vali],Reduced_Output_test[:,0], label='True Output1')
plt.plot(samples_i[N_train:-N_vali],y_pred, label='Predicted Output1')
plt.fill(np.concatenate([samples_i[N_train:-N_vali], samples_i[N_train:-N_vali][::-1]]),np.concatenate([y_pred - 3 * sigma, (y_pred + 3 * sigma)[::-1]]), alpha=.3, fc='b', ec='None')
plt.xlabel('$Sample Points$')
plt.ylabel('$Output$')
plt.legend(loc='upper right')
    
print("training error" + str(trainMSE_1))
print("testing error" + str(testMSE_1))

#%% save the model to disk
# filename = 'surrogate_model1.sav'
# pickle.dump(gp_GPML1, open(filename, 'wb'))

#%% Surrogate 2 for Output2
L0=10
L_bounds=(1e-5,1e5)

trainMSE_2=[]  
testMSE_2=[] 

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

gp_GPML2 = GaussianProcessRegressor(kernel=kernel_prod, alpha=1e-8,n_restarts_optimizer=0,normalize_y=False,random_state=123)
gp_GPML2.fit(Reduced_Input_train, Reduced_Output_train[:,1])
     
y_output = gp_GPML2.predict(Reduced_Input_train, return_std=False)
trainMSE_2.append(np.sqrt(mean_squared_error(Reduced_Output_train[:,1], y_output, multioutput='uniform_average')))

y_pred, sigma = gp_GPML2.predict(Reduced_Input_test, return_std=True)
testMSE_2.append(np.sqrt(mean_squared_error(Reduced_Output_test[:,1], y_pred, multioutput='uniform_average')))

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(samples_i[:N_train],Reduced_Output_train[:,1], label='True Output2')
plt.plot(samples_i[:N_train],y_output, label='Predicted Output2')
plt.xlabel('$Sample Points$')
plt.ylabel('$Output$')
plt.legend(loc='upper right')
plt.subplot(2, 1, 2)
plt.plot(samples_i[N_train:-N_vali],Reduced_Output_test[:,1], label='True Output2')
plt.plot(samples_i[N_train:-N_vali],y_pred, label='Predicted Output2')
plt.fill(np.concatenate([samples_i[N_train:-N_vali], samples_i[N_train:-N_vali][::-1]]),np.concatenate([y_pred - 3 * sigma, (y_pred + 3 * sigma)[::-1]]), alpha=.3, fc='b', ec='None')
plt.xlabel('$Sample Points$')
plt.ylabel('$Output$')
plt.legend(loc='upper right')
# plt.xlim(1200,1220)    
print("training error" + str(trainMSE_2))
print("testing error" + str(testMSE_2))

#%% save the model to disk
# filename = 'surrogate_model2.sav'
# pickle.dump(gp_GPML2, open(filename, 'wb'))

#%% Surrogate 3 for Output3
L0=10
L_bounds=(1e-5,1e5)

trainMSE_3=[]  
testMSE_3=[] 

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

gp_GPML3 = GaussianProcessRegressor(kernel=kernel_prod, alpha=1e-8,n_restarts_optimizer=0,normalize_y=False,random_state=612)
gp_GPML3.fit(Reduced_Input_train, Reduced_Output_train[:,2])
     
y_output = gp_GPML3.predict(Reduced_Input_train, return_std=False)
trainMSE_3.append(np.sqrt(mean_squared_error(Reduced_Output_train[:,2], y_output, multioutput='uniform_average')))

y_pred, sigma = gp_GPML3.predict(Reduced_Input_test, return_std=True)
testMSE_3.append(np.sqrt(mean_squared_error(Reduced_Output_test[:,2], y_pred, multioutput='uniform_average')))

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(samples_i[:N_train],Reduced_Output_train[:,2], label='True Output3')
plt.plot(samples_i[:N_train],y_output, label='Predicted Output3')
plt.xlabel('$Sample Points$')
plt.ylabel('$Output$')
# plt.xlim(1800,2000)
plt.legend(loc='upper right')
plt.subplot(2, 1, 2)
plt.plot(samples_i[N_train:-N_vali],Reduced_Output_test[:,2], label='True Output3')
plt.plot(samples_i[N_train:-N_vali],y_pred, label='Predicted Output3')
plt.fill(np.concatenate([samples_i[N_train:-N_vali], samples_i[N_train:-N_vali][::-1]]),np.concatenate([y_pred - 3 * sigma, (y_pred + 3 * sigma)[::-1]]), alpha=.3, fc='b', ec='None')
plt.xlabel('$Sample Points$')
plt.ylabel('$Output$')
plt.legend(loc='upper right')
    
print("training error" + str(trainMSE_3))
print("testing error" + str(testMSE_3))

#%% save the model to disk
# filename = 'surrogate_model3.sav'
# pickle.dump(gp_GPML3, open(filename, 'wb'))

#%% Surrogate 4 for Output4
L0=10
L_bounds=(1e-5,1e5)

trainMSE_4=[]  
testMSE_4=[] 

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

gp_GPML4 = GaussianProcessRegressor(kernel=kernel_prod, alpha=1e-8,n_restarts_optimizer=0,normalize_y=False,random_state=612)
gp_GPML4.fit(Reduced_Input_train, Reduced_Output_train[:,3])
     
y_output = gp_GPML4.predict(Reduced_Input_train, return_std=False)
trainMSE_4.append(np.sqrt(mean_squared_error(Reduced_Output_train[:,3], y_output, multioutput='uniform_average')))

y_pred, sigma = gp_GPML4.predict(Reduced_Input_test, return_std=True)
testMSE_4.append(np.sqrt(mean_squared_error(Reduced_Output_test[:,3], y_pred, multioutput='uniform_average')))

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(samples_i[:N_train],Reduced_Output_train[:,3], label='True Output4')
plt.plot(samples_i[:N_train],y_output, label='Predicted Output4')
plt.xlabel('$Sample Points$')
plt.ylabel('$Output$')
# plt.xlim(1800,2000)
plt.legend(loc='upper right')
plt.subplot(2, 1, 2)
plt.plot(samples_i[N_train:-N_vali],Reduced_Output_test[:,3], label='True Output4')
plt.plot(samples_i[N_train:-N_vali],y_pred, label='Predicted Output4')
plt.fill(np.concatenate([samples_i[N_train:-N_vali], samples_i[N_train:-N_vali][::-1]]),np.concatenate([y_pred - 3 * sigma, (y_pred + 3 * sigma)[::-1]]), alpha=.3, fc='b', ec='None')
plt.xlabel('$Sample Points$')
plt.ylabel('$Output$')
plt.legend(loc='upper right')
    
print("training error" + str(trainMSE_4))
print("testing error" + str(testMSE_4))

#%% save the model to disk
# filename = 'surrogate_model4.sav'
# pickle.dump(gp_GPML4, open(filename, 'wb'))


#%% Surrogate 5 for Output5
L0=10
L_bounds=(1e-5,1e5)

trainMSE_5=[]  
testMSE_5=[] 

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

gp_GPML5 = GaussianProcessRegressor(kernel=kernel_prod, alpha=1e-8,n_restarts_optimizer=0,normalize_y=False,random_state=123)
gp_GPML5.fit(Reduced_Input_train, Reduced_Output_train[:,4])
     
y_output = gp_GPML5.predict(Reduced_Input_train, return_std=False)
trainMSE_5.append(np.sqrt(mean_squared_error(Reduced_Output_train[:,4], y_output, multioutput='uniform_average')))

y_pred, sigma = gp_GPML5.predict(Reduced_Input_test, return_std=True)
testMSE_5.append(np.sqrt(mean_squared_error(Reduced_Output_test[:,4], y_pred, multioutput='uniform_average')))

plt.figure()
plt.subplot(2, 1, 1)
plt.plot(samples_i[:N_train],Reduced_Output_train[:,4], label='True Output5')
plt.plot(samples_i[:N_train],y_output, label='Predicted Output5')
plt.xlabel('$Sample Points$')
plt.ylabel('$Output$')
# plt.xlim(1800,2000)
plt.legend(loc='upper right')
plt.subplot(2, 1, 2)
plt.plot(samples_i[N_train:-N_vali],Reduced_Output_test[:,4], label='True Output5')
plt.plot(samples_i[N_train:-N_vali],y_pred, label='Predicted Output5')
plt.fill(np.concatenate([samples_i[N_train:-N_vali], samples_i[N_train:-N_vali][::-1]]),np.concatenate([y_pred - 3 * sigma, (y_pred + 3 * sigma)[::-1]]), alpha=.3, fc='b', ec='None')
plt.xlabel('$Sample Points$')
plt.ylabel('$Output$')
plt.legend(loc='upper right')
    
print("training error" + str(trainMSE_5))
print("testing error" + str(testMSE_5))

#%% save the model to disk
# filename = 'surrogate_model5.sav'
# pickle.dump(gp_GPML5, open(filename, 'wb'))


#%% Save all files in one variable
GP_list = [encoder1,decoder1,encoder2,decoder2,gp_GPML1,gp_GPML2,gp_GPML3,gp_GPML4,gp_GPML5]
file_name = "GP_Model_Global200.pkl"

pickle.dump(GP_list, open(file_name, "wb"))


#%% Load the surrogate models
 
GP_Model = pickle.load(open("GP_Model_Global200.pkl", "rb"))
encoder1 = GP_Model[0]
decoder1 = GP_Model[1]
encoder2 = GP_Model[2]
decoder2 = GP_Model[3]
Model1 = GP_Model[4]
Model2 = GP_Model[5]
Model3 = GP_Model[6]
Model4 = GP_Model[7]
Model5 = GP_Model[8]


#%% Validation
Reduced_Input_vali_0 = input_vali[:,0:nDOF] @ encoder1
Reduced_Input_vali = np.append(Reduced_Input_vali_0,np.reshape(input_vali[:,nDOF],(len(Reduced_Input_vali_0[:,1]),1)),axis=1)

Reduced_Output1 = Model1.predict(Reduced_Input_vali,return_std=False).reshape(N_vali,1)
Reduced_Output2 = Model2.predict(Reduced_Input_vali,return_std=False).reshape(N_vali,1)
Reduced_Output3 = Model3.predict(Reduced_Input_vali,return_std=False).reshape(N_vali,1)
Reduced_Output4 = Model4.predict(Reduced_Input_vali,return_std=False).reshape(N_vali,1)
Reduced_Output5 = Model5.predict(Reduced_Input_vali,return_std=False).reshape(N_vali,1)

Output_New = np.concatenate((Reduced_Output1,Reduced_Output2,Reduced_Output3,Reduced_Output4,Reduced_Output5),axis=1)
Output_New_decoded = Output_New @ decoder2