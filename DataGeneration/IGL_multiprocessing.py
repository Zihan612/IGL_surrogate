# -*- coding: utf-8 -*-
"""
Created on Mon Jan 18 16:39:10 2021

@author: MDT
"""

import sys
import os
import time
import numpy as np
from IGL import IGL


index = int(sys.argv[-3])
workPath = sys.argv[-2]
scriptPath = sys.argv[-1]
sys.path.append(scriptPath) 
cpuCount = int(sys.argv[-4])


# index = 0
# workPath = 'C:\ZIHAN\GreenupCrackEstimation\Multiprocessing\Files\\'
# workPath = workPath + "\\S" +str(index) + "\\"
# scriptPath = 'C:\ZIHAN\GreenupCrackEstimation\Multiprocessing\\'
# sys.path.append(scriptPath) 
# cpuCount = 4


#%%
IGLInfo = {}
globInfo = {}
locInfo = {}


IGLInfo['globType'] = 'Abaqus'# 'Abaqus' or 'AbaqusSubStructure' 
IGLInfo['locType'] = 'Abaqus' # 'Abaqus', 'AbaqusSubStructure', or 'GaussianProcess'
IGLInfo['algorithm'] = 'FixedPointAitken' # 'FixedPoint', 'FixedPointAitken', or 'ConjugateGradient
IGLInfo['tol'] = 10**(-5)
IGLInfo['maxSteps'] = 20
IGLInfo['SIFConvergence'] = True
IGLInfo['abaqusVersion'] = 'abq2021'
IGLInfo['CPU'] = cpuCount


globInfo['saveResults'] = True #
globInfo['instanceName'] = 'Part-1-failed-1'#'Part-1-1'
globInfo['partNames'] = 'Part-1-failed'#'Part-1'
globInfo['boundary'] = 'LOCBOUNDARY' # Assembly set names. Must be all capitalized. Local boundary in global cae
globInfo['globLocDomain'] = 'LOCDOMAIN'
globInfo['stepName'] = 'Step-1'
globInfo['CAEName'] = 'GreenupGlobal'+str(index)+'.cae'
globInfo['modelName'] = 'GlobModel-'+str(index)
globInfo['SSSavedFlag'] = True
globInfo['cpuCount'] = cpuCount



# for local type = Abaqus
# Assumed that local solid domain has 'Solid' in part name
locInfo['saveResults'] = True #
locInfo['instanceName'] = 'Part-1-failed-1'#'Part-1-1'
locInfo['partNames'] = 'Part-1-failed'#'Part-1'
locInfo['crackPartName'] = 'Crack'
locInfo['boundary'] = 'LOCBOUNDARY' # Assembly set names. Must be all capitalized. Local boundary in global cae
locInfo['globLocDomain'] = 'LOCDOMAIN'
locInfo['stepName'] = 'Step-1'
locInfo['CAEName'] = 'GreenupLocal'+str(index)+'.cae'
locInfo['abaqusVersion'] = 'abq2021'
locInfo['cpuCount'] = cpuCount
locInfo['thickness'] = 0.75
locInfo['normal'] = np.array([0.,-1.,0.])
locInfo['propDir'] = np.array([-1.,0.,0.]) # Need to figure this out
locInfo['translation'] = np.array([-284.,-212.5,-723.]) # x,y,z coordinates of center of plate at which crack starts

Job = {}
Job['identifier'] = 0
Job['index'] = index
Job['globInfo'] = {}
Job['locInfo'] = {}
Job['IGLInfo'] = {}
Job['scriptPath']=scriptPath
Job['locPath'] = workPath
Job['saveResults'] = True

try:
    os.mkdir(Job['locPath'])
except:
    pass
for key,value in IGLInfo.items():
    Job['IGLInfo'][key] = value
for key,value in globInfo.items():
    Job['globInfo'][key] = value
for key,value in locInfo.items():
    Job['locInfo'][key] = value

Job['locInfo']['baseCAEName'] = locInfo['CAEName']
Job['locInfo']['CAEName'] = locInfo['CAEName']

t2 = time.time()
IGLInstance = IGL(Job)
t3 = time.time()
print('IGL Initiation time is '+str(t3-t2))
IGLInstance.run()
t4 = time.time()
