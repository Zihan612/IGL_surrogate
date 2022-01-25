# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 10:21:37 2021

@author: MDT
"""

from odbAccess import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *
import pickle
import numpy as np
import csv
import os.path
import sys
# import pandas as pd

index = int(sys.argv[-3])
# index = 0
workPath = sys.argv[-2]
# workPath = 'C:\ZIHAN\GreenupCrackEstimation\Multiprocessing\Files\S0'
scriptPath = sys.argv[-1]
# scriptPath = 'C:\ZIHAN\GreenupCrackEstimation\Multiprocessing'

os.chdir(workPath)

## CHANGE odbName to the name of your ODB file

odbName='GreenupLocal'+str(index)+'_postprocess.odb'
# Opening the Odb File
odb = openOdb(workPath+'\\'+ odbName)


## Create the CSV file based on "odbName" ~~ odbNameStrain.csv. If the file already exists, overwrite it
fileName=scriptPath + '\\Files\OutputData\SIF\S{}.csv'.format(index)
		
#Get the desired SIF values

allOutputs = odb.steps['Step-1'].historyRegions['ElementSet  ALL ELEMENTS'].historyOutputs
for ele in range(1,14):
    SIFValues=np.zeros(3) #creates an empty list that you can dynamically append SIF values to
    for Ks in range(1,4):
        SIFName = 'K{} at CRACK-1_XFEM_{}_Contour_4'.format(Ks,ele)
        SIF = allOutputs[SIFName].data[0]
        SIFValues[Ks-1] = SIF[-1]#append the strain for this step and set to the dynamic list
    with open(fileName, 'ab') as csvfile:
    		fileWriter=csv.writer(csvfile)
    		fileWriter.writerow(SIFValues)

# Closing Odb
odb.close()  #this closes the ODB file. Not necessary, but if you want it closed, remove the comment marker