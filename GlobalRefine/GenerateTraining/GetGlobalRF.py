# -*- coding: utf-8 -*-
"""
Created on Thu Sep 10 08:43:08 2020

@author: RDCHLTBF
"""

from odbAccess import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *
import numpy as np
import sys
import pickle

# Import variables from driver
initialRF = False
if sys.argv[-6] == 'True':
    initialRF = True
elif sys.argv[-6] != 'False':
    sys.exit('Error: wrong input for initialRF')
fileOutputName = sys.argv[-7]
setLocal = sys.argv[-5]
setNames = sys.argv[-4]
stepName = sys.argv[-3]
odbName = str(sys.argv[-2])
workPath = str(sys.argv[-1])

fPrint = open(workPath+'\\stderr.txt','a') 
sys.stdout=fPrint

odb = openOdb(str(workPath+'\\'+ odbName), readOnly=FALSE) # Open the Odb File
Frame=odb.steps[stepName].frames[-1]  # There is only one step

# Get all nodes in setNames for global
AbRF = {}
AbNodes = {}
AbLabels = []
AbNodesLabel = {}
AbConnect = {}
AbConnectInd = {}
AbElems = {}

if initialRF == True:
    
    #print(odb.rootAssembly.elementSets[iSet])  # List of all nodes along local boundary
    AbElems[setNames] = []
    AbElems[setNames] = odb.rootAssembly.elementSets[setLocal].elements[0]
    AbConnect[setNames] = []
    AbConnectInd[setNames] = []
    for j in odb.rootAssembly.elementSets[setLocal].elements[0]:
        AbConnect[setNames] += [[]]
        AbConnectInd[setNames] += [[]]
    
    AbNodes[setNames] = [] # odb.rootAssembly.nodeSets[setNames].nodes[0]  # List of all nodes along local boundary
    AbLabels = np.zeros(len(odb.rootAssembly.nodeSets[setNames].nodes[0])) # Standardize the order of the array by node labels
    for x, n in enumerate(odb.rootAssembly.nodeSets[setNames].nodes[0]):
        AbLabels[x] = n.label
    mapping = np.argsort(AbLabels,axis=0)
    for i,iMap in enumerate(mapping):
        AbNodes[setNames] += [odb.rootAssembly.nodeSets[setNames].nodes[0][iMap]]
    
    # This is to help debug, not strictly necessary and may be removed for speed
    for j,jNodes in enumerate(AbNodes[setNames]):
        for k,kEls in enumerate(AbElems[setNames]):
            if jNodes.label in kEls.connectivity:
                AbConnect[setNames][k] += [jNodes.label]
                AbConnectInd[setNames][k] += [j]
    AbRF[setNames] = []
    
    #filename = 'RFElems'
    #outfile = open(filename,'wb')
    #pickle.dump(AbElems,outfile)
    #outfile.close()
    filename = 'RFConnect'
    outfile = open(workPath+'\\'+filename,'wb')
    pickle.dump(AbConnect,outfile)
    outfile.close()
    filename = 'RFConnectInd'
    outfile = open(workPath+'\\'+filename,'wb')
    pickle.dump(AbConnectInd,outfile)
    outfile.close()
    print('Finished finding elements and nodes for global reaction')
else:
    AbElems[setNames] = []
    AbElems[setNames] = odb.rootAssembly.elementSets[setLocal].elements[0]
    filename = 'RFConnect'
    infile = open(workPath+'\\'+filename,'rb')
    AbConnect = pickle.load(infile)
    infile.close()
    filename = 'RFConnectInd'
    infile = open(workPath+'\\'+filename,'rb')
    AbConnectInd = pickle.load(infile)
    infile.close()
    
    AbNodes[setNames] = [] # odb.rootAssembly.nodeSets[setNames].nodes[0]  # List of all nodes along local boundary
    AbLabels = np.zeros(len(odb.rootAssembly.nodeSets[setNames].nodes[0])) # Standardize the order of the array by node labels
    for x, n in enumerate(odb.rootAssembly.nodeSets[setNames].nodes[0]):
        AbLabels[x] = n.label
    mapping = np.argsort(AbLabels,axis=0)
    for i,iMap in enumerate(mapping):
        AbNodes[setNames] += [odb.rootAssembly.nodeSets[setNames].nodes[0][iMap]]
    
    AbRF[setNames] = []
    

ResultField1 = Frame.fieldOutputs['NFORC1']
ResultField2 = Frame.fieldOutputs['NFORC2']
ResultField3 = Frame.fieldOutputs['NFORC3']
ResultField4 = Frame.fieldOutputs['NFORC4']
ResultField5 = Frame.fieldOutputs['NFORC5']
ResultField6 = Frame.fieldOutputs['NFORC6']

for key,values in AbElems.items():
    iElSet = values
    AbRF[key] = np.zeros((len(AbNodes[key]),6))
    for j,jEl in enumerate(AbConnect[key]):
        if not jEl:
            continue
        else:
            
            jAbEl = iElSet[j]
            jNFORC1 = ResultField1.getSubset(region=jAbEl).values
            jNFORC2 = ResultField2.getSubset(region=jAbEl).values
            jNFORC3 = ResultField3.getSubset(region=jAbEl).values
            jNFORC4 = ResultField4.getSubset(region=jAbEl).values
            jNFORC5 = ResultField5.getSubset(region=jAbEl).values
            jNFORC6 = ResultField6.getSubset(region=jAbEl).values
            
            for k,kNFORC1 in enumerate(jNFORC1):
                if kNFORC1.nodeLabel in AbConnect[key][j]:
                    kInd = AbConnect[key][j].index(kNFORC1.nodeLabel)
                    AbRF[key][AbConnectInd[key][j][kInd],0] += kNFORC1.dataDouble
                    AbRF[key][AbConnectInd[key][j][kInd],1] += jNFORC2[k].dataDouble
                    AbRF[key][AbConnectInd[key][j][kInd],2] += jNFORC3[k].dataDouble
                    AbRF[key][AbConnectInd[key][j][kInd],3] += jNFORC4[k].dataDouble
                    AbRF[key][AbConnectInd[key][j][kInd],4] += jNFORC5[k].dataDouble
                    AbRF[key][AbConnectInd[key][j][kInd],5] += jNFORC6[k].dataDouble
odb.close()
# use npy
np.save(workPath+'\\'+fileOutputName+'RF',AbRF[setNames])