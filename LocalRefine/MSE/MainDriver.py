# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 14:55:41 2021

@author: MDT
"""

import sys
import os

import shutil
import glob

from subprocess import check_call
import subprocess

import pickle


#%% Set the scripPath and workPath

scriptPath = os.path.dirname(os.path.abspath(__file__)) + "\\"
os.chdir(scriptPath)
sys.path.append(scriptPath)

workPath= scriptPath + "Files\\"
if not os.path.exists(workPath):
    os.makedirs(workPath)

#%% Pre-defined function

# use os.chdir( ) first to define the path
def local_new(name,path): # Create a new cae file to modify
    s = scriptPath + "\\GreenupLocal.cae"
    shutil.copy(s, path + "\\GreenupLocal" + name +".cae")
    return None

# use os.chdir( ) first to define the path
def clear_files(option): # Clear the residual files 1: keep odb; 0: delete odb

    fileList = glob.glob('*.*')
    files_keep =  ['ExtractDispl_CoordVals.py','Reaction.sav','TrainingSpace.sav','Input_train_Global200.sav','Output_train_Global200.sav', 'UpdatingSurrogate.py', 'AddedPoints.sav', 'InputDecoder.sav','InputEncoder.sav','OutputEncoder.sav','OutputDecoder.sav','MC_space.sav','MC_denorm.sav', 'UpdatingDriver.py', 'IGL_multiprocessing.py', 'AbaqusEditMdb.py', 'localFEM.py', 'GetOdbResults.py', 'GetOdbSIFResults.py','MainDriver.py', 'ShellAbaqus.py', 'SubmitAbaqusJob.py','GreenupLocal.cae']
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

def FEM_run(i,script_directory,working_directory,cpuCount):
      
    
    os.chdir(script_directory)
    try:
        check_call([sys.executable,"UpdatingDriver.py",str(i),str(script_directory)], stdin=None, stdout=None, stderr=None, shell=True, timeout =50000.0)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        return 1
    
    working_directory = working_directory + "\\Model" +str(i+1) +"\\"
    os.chdir(working_directory)
    PointUpdating = pickle.load(open('AddedPoints.sav', 'rb'))
    
    for N_point in range(len(PointUpdating[:,0])):
        localPath = working_directory + "\\S" +str(N_point) +"\\"
        if not os.path.exists(localPath):
            os.makedirs(localPath) 
        local_new(str(N_point),localPath)
    
        os.chdir(script_directory)
        try:
            check_call("abq2021 cae noGUI=localFEM.py -- " +str(N_point)+ " " +str(cpuCount)+ " "+str(localPath)+ " " +str(PointUpdating[N_point,-1]), stdin=None, stdout=None, stderr=None, shell=True, timeout = 500.0)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
            return 1
     
    os.chdir(script_directory)
    try:
        check_call([sys.executable,"IGL_multiprocessing.py",str(i),str(cpuCount),str(working_directory),str(script_directory),str(len(PointUpdating[:,0]))], stdin=None, stdout=None, stderr=None, shell=True, timeout =50000.0)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        return 1
    
    os.chdir(script_directory)
    try:
        check_call([sys.executable,"UpdatingSurrogate.py",str(i),str(working_directory),str(script_directory)], stdin=None, stdout=None, stderr=None, shell=True, timeout =50000.0)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
        return 1
    
    os.chdir(script_directory)
    clear_files(1)
    # os.chdir(working_directory)
    # clear_files(0)

os.chdir(scriptPath)
clear_files(1)

#%% Single Run

for index in range(100):
    i,script_directory,working_directory,cpuCount = [index,scriptPath,workPath,4]
    FEM_run(i,script_directory,working_directory,cpuCount)

#%% Test
# index=5
# i,script_directory,working_directory,cpuCount = [index,scriptPath,workPath,4]
# FEM_run(i,script_directory,working_directory,cpuCount)

