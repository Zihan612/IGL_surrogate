# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 14:09:26 2020

@author: RDCHLTBF
"""

from odbAccess import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *
import numpy as np
import sys

    
# Import variables from driver
fileResultNames = str.strip(sys.argv[-6])
fileOutputName = str.strip(sys.argv[-5])
setNames = str(sys.argv[-4])
stepName = str(sys.argv[-3])
odbName = str(sys.argv[-2])
workPath = str(sys.argv[-1])

fPrint = open(workPath+'\\stderr.txt','a') 
sys.stdout=fPrint
        
resultNames = []
with open(workPath+'\\'+fileResultNames) as f2:
    for cnt, line in enumerate(f2):
        resultNames += [str.strip(line)]

odb = openOdb(str(workPath+'\\'+ odbName), readOnly=FALSE) # Open the Odb File
print(str(workPath+'\\'+ odbName))
print(stepName)
print(odb)
print(odb.steps)
print(odb.steps[stepName])
print(odb.steps[stepName].frames)
Frame=odb.steps[stepName].frames[-1]  # There is only one step

# Get all nodes in setNames for global
AbNodes = {}
AbResults = {}

AbNodes[setNames] = [] # odb.rootAssembly.nodeSets[setNames].nodes[0]  # List of all nodes along local boundary
AbLabels = np.zeros(len(odb.rootAssembly.nodeSets[setNames].nodes[0])) # Standardize the order of the array by node labels
for x, n in enumerate(odb.rootAssembly.nodeSets[setNames].nodes[0]):
    AbLabels[x] = n.label
mapping = np.argsort(AbLabels,axis=0)
for i,iMap in enumerate(mapping):
    AbNodes[setNames] += [odb.rootAssembly.nodeSets[setNames].nodes[0][iMap]]

AbResults[setNames] = {}
for j,jResult in enumerate(resultNames):
    AbResults[setNames][jResult] = []
for iR in resultNames:
    if iR == 'Coord':
                
        for key,values in AbNodes.items():
            iSet = values
            AbResults[key][iR] = np.zeros((len(iSet),3))

            for j,jNode in enumerate(iSet):
                AbResults[key][iR][j,:]   = jNode.coordinates
                
    else:
        ResultField = Frame.fieldOutputs[iR]
        
        for key,values in AbNodes.items():
            iSet = values
            AbResults[key][iR] = np.zeros((len(iSet),3))
            for j,jNode in enumerate(iSet):
                AbResults[key][iR][j,:] = ResultField.getSubset(region=jNode).values[0].dataDouble
odb.close()

# use npy
for iR in resultNames:
    np.save(workPath+'\\'+iR+fileOutputName,AbResults[setNames][iR])