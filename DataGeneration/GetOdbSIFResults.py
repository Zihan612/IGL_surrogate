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
fileOutputName = sys.argv[-4]
stepName = str(sys.argv[-3])
odbName = str(sys.argv[-2])
workPath = str(sys.argv[-1])
fPrint = open(workPath+'\\stderr.txt','a') 
sys.stdout=fPrint


#stepName = 'Step-1'
#workPath = 'C:\\Temp\\'
#odbName = 'Local_Toy_postprocess.odb'
#odbName = 'Local_Toy_job6.odb'
#odbName = 'Toy_reference.odb'
#odbName = 'GreenupGlobal_reference.odb'
#odbName = 'GreenupLocalQuarter_job4.odb'



odb = openOdb(str(workPath+'\\'+ odbName), readOnly=FALSE) # Open the Odb File
print(str(workPath+'\\'+ odbName))
print(stepName)
print(odb)
print(odb.steps)
print(odb.steps[stepName])
print(odb.steps[stepName].frames)
Frame=odb.steps[stepName].frames[-1]  # There is only one step

# Get all nodes in setNames for global
AbResults = []
for key,values in odb.steps[stepName].historyRegions.items():
    for key2,values2 in values.historyOutputs.items():
        if 'K1' in key2 and 'Contour_4' in key2:
            AbResults += [values2.data[-1][1]]
#print(Frame.fieldOutputs)
AbResults = np.array(AbResults)
print(AbResults)
odb.close()

np.save(workPath+'\\SIF'+fileOutputName,AbResults)
print(AbResults)