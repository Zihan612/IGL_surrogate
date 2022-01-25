# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 16:39:10 2021

@author: MDT
"""

import sys
import os
import numpy as np
from ShellAbaqus import ShellAbaqus

import pickle


#%%
index = int(sys.argv[-5])

working_directory = sys.argv[-3]

scriptPath = sys.argv[-2]

sys.path.append(scriptPath)
 
cpuCount = int(sys.argv[-4])

N_runs = int(sys.argv[-1])


# index = 0
# working_directory = 'C:\ZIHAN\GreenupCrackEstimation\A1201\LocalRefine\MSE\Files\Model1\\'
# scriptPath = 'C:\ZIHAN\GreenupCrackEstimation\A1201\LocalRefine\MSE\\'
# sys.path.append(scriptPath)
# cpuCount = 4
# N_runs = 2


#%%
os.chdir(working_directory)
InputSpace0 = pickle.load(open('AddedPoints.sav', 'rb'))
InputSpace = InputSpace0[:,0:4]

os.chdir(scriptPath)
InputDecoder = pickle.load(open('InputDecoder.sav', 'rb'))
InputReconstruct = InputSpace @ InputDecoder

OutputEncoder = pickle.load(open('OutputEncoder.sav', 'rb'))

OutputReaction = []



#%%
for runs in range(N_runs):
    abaqusVersion = 'abq2021'
    workPath = working_directory + "\\S" +str(runs) +"\\"

    locInfo = {}
    locInfo['saveResults'] = True #
    locInfo['instanceName'] = 'Part-1-failed-1'#'Part-1-1'
    locInfo['partNames'] = 'Part-1-failed'#'Part-1'
    locInfo['crackPartName'] = 'Crack'
    locInfo['boundary'] = 'LOCBOUNDARY' # Assembly set names. Must be all capitalized. Local boundary in global cae
    locInfo['globLocDomain'] = 'LOCDOMAIN'
    locInfo['stepName'] = 'Step-1'
    locInfo['CAEName'] = 'GreenupLocal'+str(runs)+'.cae'
    locInfo['abaqusVersion'] = 'abq2021'
    locInfo['cpuCount'] = cpuCount
    locInfo['thickness'] = 0.75
    locInfo['normal'] = np.array([0.,-1.,0.])
    locInfo['propDir'] = np.array([-1.,0.,0.]) # Need to figure this out
    locInfo['translation'] = np.array([-284.,-212.5,-723.]) # x,y,z coordinates of center of plate at which crack starts
    locInfo['baseCAEName'] = locInfo['CAEName']
    locInfo['CAEName'] = locInfo['CAEName']
    locInfo['inpName'] = locInfo['CAEName'][:-4]+'.inp'
                
    
    workDir = os.listdir(workPath)
    for item in workDir: # Remove all lock files in directory
        if item.endswith(".lck"):
            os.remove(os.path.join(workPath, item))
    workDir = os.listdir(workPath)
    for item in workDir: # Remove all job files in directory
        if 'job' in item: 
            os.remove(os.path.join(workPath, item))
    
    Loc = ShellAbaqus(abaqusVersion,workPath,scriptPath,locInfo)
    
    
    # FixedPointAitken()
    iteration = 0
    
    LocalAppliedDispl = {}
    LocalReactionTrain = {}
    
    while iteration == 0:
        
        LocalAppliedDispl_temp = np.reshape(InputReconstruct[runs,:],(120,6))
        LocalAppliedDispl[iteration] = np.array(LocalAppliedDispl_temp)
        Loc.apply_displacement(LocalAppliedDispl[iteration],iteration)
        Loc.run_job(iteration)
        LocalReactionTrain[iteration] = Loc.get_odb_RF(iteration)
        Reaction_temp = np.reshape(np.array(LocalReactionTrain[iteration]),(120*6,))
        ReducedReaction = Reaction_temp @ OutputEncoder
        OutputReaction.append(ReducedReaction)
        iteration += 1
        
        
#%%
os.chdir(working_directory)
filename = 'Reaction.sav'
pickle.dump(OutputReaction, open(filename, 'wb'))