# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 14:55:41 2021

@author: MDT
"""

import sys
import os

import shutil
import glob

import pickle

from multiprocessing import Pool
import time

from subprocess import check_call
import subprocess


start_time = time.time()

#%% Set the scripPath and workPath

scriptPath = os.path.dirname(os.path.abspath(__file__)) + "\\"
os.chdir(scriptPath)
sys.path.append(scriptPath)

workPath= scriptPath + "Files\\"
if not os.path.exists(workPath):
    os.makedirs(workPath)
    
dataPath= workPath + "OutputData\\"
if not os.path.exists(dataPath):
    os.makedirs(dataPath)


#%% Pre-defined function

def local_new(name,path): # Create a new cae file to modify
    s = scriptPath + "\\GreenupLocal.cae"
    shutil.copy(s, path + "\\GreenupLocal" + name +".cae")
    return None


# use os.chdir( ) first to define the path
def clear_files(option): # Clear the residual files 1: keep odb; 0: delete odb


    fileList = glob.glob('*.*')
    files_keep =  ['InputDecoder.sav','FillingPoints300.sav','ExtractDispl_CoordVals.py','SIFExtractor.py', 'IGL_multiprocessing.py','IGL.py', 'AbaqusEditMdb.py', 'localFEM.py', 'globalFEM.py', 'GetOdbResults.py', 'GetOdbSIFResults.py','MainDriver.py', 'ShellAbaqus.py','SubmitAbaqusJob.py','GetGlobalRF.py','GreenupGlobal.cae','GreenupLocal.cae']
    if option == 1:
        odb_files = glob.glob('*.odb')
        # cae_files = glob.glob('*.cae')
        files_keep = files_keep + odb_files
    else:
        files_keep = files_keep
    
    for ele in files_keep:
        try:
            fileList.remove(ele)
        except (ValueError) as e:
            pass
    
    for item in fileList:
        try:
            os.remove(item)
        except OSError:
            pass
        
    return None


#%% For single run

def FEM_run(i,working_directory,script_directory, crack_length, cpuCount):
    
    
################################################################################
# CRACK FEM
################################################################################    
    working_directory= working_directory + "\\S" +str(i) +"\\"
    if not os.path.exists(working_directory):
        os.makedirs(working_directory)
    
    local_new(str(i),working_directory)

    os.chdir(script_directory)
    try:
        check_call("abq2021 cae noGUI=localFEM.py -- " +str(cpuCount)+ " " +str(i)+ " " +str(working_directory)+ " " +str(crack_length), stdin=None, stdout=None, stderr=None, shell=True, timeout = 500.0)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        return 1   

        
    os.chdir(script_directory)
    try:
        check_call([sys.executable,"IGL_multiprocessing.py",str(cpuCount),str(i),str(working_directory),str(script_directory)], stdin=None, stdout=None, stderr=None, shell=True, timeout =5000.0)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        return 1
    
    os.chdir(script_directory)
    clear_files(1)
    os.chdir(working_directory)
    clear_files(0)

os.chdir(scriptPath)
clear_files(1)

#%% Single Run
#### Full length: 0~8; Half length: 0~3.75####
# h_up and h_down must be multiples of 24
# n_samples = 200
# distribution = chaospy.Iid(chaospy.Uniform(0, 1), 4)
# samples = distribution.sample(n_samples,rule="sobol")
# Input4 = samples[3,:]*3.75


# index=0
# i,working_directory,script_directory,crack_length,cpuCount = [index,workPath,scriptPath,Input4[index,],2]
# FEM_run(i,working_directory,script_directory,crack_length, cpuCount)

#%% Monte Carlo

# MissingList = [3,8]

InputSpace = pickle.load(open('FillingPoints300.sav', 'rb'))
Input4 = InputSpace[:,-1]
cpuCount = 2 # How many cpus you want to use for a single process
  
CPU = 2 # how many processes run in parallel
#%%

if __name__ == "__main__":
    pool = Pool(processes=CPU)
    
    time.sleep(3)
    
    results = pool.starmap(FEM_run,[(i,workPath,scriptPath,Input4[i,],cpuCount) for i in range(300)],)
    # results = pool.starmap(FEM_run,[(i,workPath,scriptPath,Input4[i,],cpuCount) for i in MissingList],)
    
    pool.close()
    
    print("--- Execution time: %s seconds ---" % (time.time() - start_time))


#%%
# SIFPath = "C:\ZIHAN\GreenupCrackEstimation\SequentialSampling\Files\OutputData\\"
# N_total = 300

# indexList = []

# for index in range(0,N_total):
#     fileName = SIFPath + "S" + str(index) + "_Reaction.csv"
#     if os.path.exists(fileName):
#         indexList += [index]

# MissingList = []
# for index in range(0,N_total):
#     if index not in indexList:
#         MissingList.append(index)



