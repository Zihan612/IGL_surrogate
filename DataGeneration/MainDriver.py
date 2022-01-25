# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 14:55:41 2021

@author: MDT
"""

import sys
import os

import shutil
import glob

import chaospy

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
    
displPath= dataPath + "LocalAppliedDispl\\"
if not os.path.exists(displPath):
    os.makedirs(displPath)
    
reactionPath= dataPath + "LocalReaction\\"
if not os.path.exists(reactionPath):
    os.makedirs(reactionPath)
    
SIFPath= dataPath + "SIF\\"
if not os.path.exists(SIFPath):
    os.makedirs(SIFPath)


#%% Pre-defined function

# use os.chdir( ) first to define the path
def global_new(name,path): # Create a new cae file to modify
    s = scriptPath + "\\GreenupGlobal.cae"
    shutil.copy(s, path + "\\GreenupGlobal" + name +".cae")
    return None

def local_new(name,path): # Create a new cae file to modify
    s = scriptPath + "\\GreenupLocal.cae"
    shutil.copy(s, path + "\\GreenupLocal" + name +".cae")
    return None

def aux_new(name,path): # Create a new cae file to modify
    s = scriptPath + "\\GreenupGlobal_aux.cae"
    shutil.copy(s, path + "\\GreenupGlobal" + name +"_aux.cae")
    return None

# use os.chdir( ) first to define the path
def clear_files(option): # Clear the residual files 1: keep odb; 0: delete odb


    fileList = glob.glob('*.*')
    files_keep =  ['ExtractDispl_CoordVals.py','SIFExtractor.py', 'IGL_multiprocessing.py','IGL.py', 'AbaqusEditMdb.py', 'localFEM.py', 'globalFEM.py', 'GetOdbResults.py', 'GetOdbSIFResults.py','MainDriver.py', 'ShellAbaqus.py','SubmitAbaqusJob.py','GetGlobalRF.py','GreenupGlobal.cae','GreenupLocal.cae']
    if option == 1:
        odb_files = glob.glob('*_postprocess.odb')
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

def FEM_run(i,working_directory,script_directory, h_up, h_down, gap_length, crack_length, cpuCount):
    
    
################################################################################
# CRACK FEM
################################################################################    
    working_directory= working_directory + "\\S" +str(i) +"\\"
    if not os.path.exists(working_directory):
        os.makedirs(working_directory)
    
    local_new(str(i),working_directory)
    global_new(str(i),working_directory)

    os.chdir(script_directory)
    try:
        check_call("abq2021 cae noGUI=localFEM.py -- " +str(cpuCount)+ " " +str(i)+ " " +str(working_directory)+ " " +str(crack_length), stdin=None, stdout=None, stderr=None, shell=True, timeout = 500.0)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        return 1   

    os.chdir(script_directory)
    try:
        check_call("abq2021 cae noGUI=globalFEM.py -- " +str(cpuCount)+ " " +str(i)+ " " +str(working_directory)+ " " +str(h_up)+ " " +str(h_down)+ " " +str(gap_length), stdin=None, stdout=None, stderr=None, shell=True, timeout = 500.0)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        return 1
        
    os.chdir(script_directory)
    try:
        check_call([sys.executable,"IGL_multiprocessing.py",str(cpuCount),str(i),str(working_directory),str(script_directory)], stdin=None, stdout=None, stderr=None, shell=True, timeout =5000.0)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        return 1
    
    os.chdir(script_directory)
    try:
        check_call("abq2021 cae noGUI=SIFExtractor.py -- " +str(i)+ " " +str(working_directory)+ " " +str(script_directory), stdin=None, stdout=None, stderr=None, shell=True, timeout = 500.0)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        return 1
    os.chdir(script_directory)
    clear_files(1)
    os.chdir(working_directory)
    clear_files(1)

os.chdir(scriptPath)
clear_files(1)

#%% Single Run
#### Full length: 0~8; Half length: 0~3.75####
# h_up and h_down must be multiples of 24

n_samples = 400
distribution = chaospy.Iid(chaospy.Uniform(0, 1), 4)
samples = distribution.sample(n_samples,rule="sobol")
Input1 = 432+samples[0,:]*288 # range:36*12-60*12
Input2 = 120+samples[1,:]*240 # range:10*12-30*12
Input3 = samples[2,:]*150 # range:10*12-30*12
Input4 = 0.5+samples[3,:]*3.5


# index=2
# i,working_directory,script_directory, h_up, h_down, gap_length, crack_length, cpuCount = [index,workPath,scriptPath,Input1[index,],Input2[index,],Input3[index,],Input4[index,],8]
# FEM_run(i,working_directory,script_directory, h_up, h_down, gap_length,crack_length, cpuCount)

#%% Monte Carlo

# MissingList = [3,8]

cpuCount = 2 # How many cpus you want to use for a single process
  
CPU = 2 # how many processes run in parallel

if __name__ == "__main__":
    pool = Pool(processes=CPU)
    
    time.sleep(3)
    
    results = pool.starmap(FEM_run,[(i,workPath,scriptPath,Input1[i,],Input2[i,],Input3[i,],Input4[i,],cpuCount) for i in range(5)],)
    # results = pool.starmap(FEM_run,[(i,workPath,scriptPath,Input1[i,],Input2[i,],Input3[i,],Input4[i,],cpuCount) for i in MissingList],)
    
    pool.close()
    
    print("--- Execution time: %s seconds ---" % (time.time() - start_time))


#%%
# SIFPath = "C:\ZIHAN\GreenupCrackEstimation\A1201\DataGeneration\Files\OutputData\SIF\\"
# N_total = 400

# indexList = []

# for index in range(0,N_total):
#     fileName = SIFPath + "S" + str(index) + ".csv"
#     if os.path.exists(fileName):
#         indexList += [index]

# MissingList = []
# for index in range(0,N_total):
#     if index not in indexList:
#         MissingList.append(index)



