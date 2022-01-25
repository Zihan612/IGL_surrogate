# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 10:19:29 2021

@author: MDT
"""
import os
import sys
import pandas as pd
import numpy as np

from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import Matern, ConstantKernel as C_ker

import pickle

#%% Passing variables

index = int(sys.argv[-3])
# index=4

workPath = sys.argv[-2]
# workPath = 'C:\\ZIHAN\\MiterGateAbauqs\\AASurrogateRefine\\AdaptiveLocalRefine\\Files\Model5\\'
# working_directory = 'C:\\ZIHAN\\MiterGateAbauqs\\AASurrogateRefine\\AdaptiveLocalRefine\\Files\\Model1\\'

scriptPath = sys.argv[-1]
# scriptPath = 'C:\\ZIHAN\\MiterGateAbauqs\\AASurrogateRefine\\AdaptiveLocalRefine\\'
# scriptPath = 'C:\\ZIHAN\\MiterGateAbauqs\\AASurrogateRefine\\AdaptiveLocalRefine\\'
sys.path.append(scriptPath)

#%% Load previous training dataset
lastModelPath = scriptPath + "\\Files\\Model" +str(index) +"\\"
os.chdir(lastModelPath)
Input0 = pickle.load(open('Input_train_Global200.sav', 'rb'))
Output0 = pickle.load(open('Output_train_Global200.sav', 'rb'))

#%% Add New Training Samples
#%% Read New Output Space
os.chdir(workPath)
OutputSpace = pickle.load(open('Reaction.sav', 'rb'))
Reduced_Output_train = np.append(Output0,OutputSpace,axis=0)
filename = 'Output_train_Global200.sav'
pickle.dump(Reduced_Output_train, open(filename, 'wb')) 

#%% Read New Input Space
InputSpace = pickle.load(open('AddedPoints.sav', 'rb'))
Reduced_Input_train = np.append(Input0,InputSpace,axis=0)
filename = 'Input_train_Global200.sav'
pickle.dump(Reduced_Input_train, open(filename, 'wb'))                               