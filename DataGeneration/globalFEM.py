# -*- coding: utf-8 -*-
"""
Created on Thu Jan 14 16:26:40 2021

@author: MDT
"""
from abaqus import *
from abaqusConstants import *
import __main__

import sys
import os

import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import connectorBehavior

cpuCount = int(sys.argv[-6])
index = int(sys.argv[-5])

workPath = sys.argv[-4]
os.chdir(workPath)

h_up = float(sys.argv[-3])
h_down = float(sys.argv[-2])

gap_length = float(sys.argv[-1])

hblock = 762.

glnorm = gap_length/hblock


namejob = 'GreenupGlobal'+str(index)
mdb=openMdb(namejob+ '.cae')

namemodel = 'GlobModel-'+str(index)
mdb.models.changeKey(fromName='Model-1', toName=namemodel)

a = mdb.models[namemodel].rootAssembly

# Gap Length
PartContact = mdb.models[namemodel].parts['ContactBlock']

e = PartContact.edges
pickedEdges = e.findAt(((-45.25, 0.0, 190.5), ))
PartContact.PartitionEdgeByParam(edges=pickedEdges, parameter=glnorm)

a.regenerate()

e1 = a.instances['ContactBlock-1'].edges
edges1 = e1.findAt(((-656.692295, -208.897428, -20), ))
region = regionToolset.Region(edges=edges1)
mdb.models[namemodel].DisplacementBC(name='QuoinBC', createStepName='Step-1', 
        region=region, u1=0.0, u2=0.0, u3=UNSET, ur1=UNSET, ur2=UNSET, 
        ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, 
        fieldName='', localCsys=None)

PartContact.generateMesh()


# Hydraulic pressures
botPoint = -777
region = a.surfaces['Upstream']

# x_up = int(h_up/12.) # Upstream hydraulic pressure
# if (x_up< 10.):
#     loadString = 'Water0'+str(x_up)
# else:
#     loadString = 'Water'+str(x_up)
    
# mdb.models[namemodel].Pressure(name=loadString, 
#     createStepName='Step-1', region=region, distributionType=HYDROSTATIC, 
#     field='', magnitude=(x_up*62.4/1000./144.), hZero=botPoint+(x_up*12.), hReference = botPoint)

# x_down = int(h_down/12.) # Downstream hydraulic pressure
# if (x_down< 10.):
#     loadString = 'Water0'+str(x_down)
# else:
#     loadString = 'Water'+str(x_down)
    
# mdb.models[namemodel].Pressure(name=loadString, 
#     createStepName='Step-1', region=region, distributionType=HYDROSTATIC, 
#     field='', magnitude=(-x_down*62.4/1000./144.), hZero=botPoint+(x_down*12.), hReference = botPoint)

x_up = int(h_up/12.) # Upstream hydraulic pressure

    
mdb.models[namemodel].Pressure(name='Upstream', 
    createStepName='Step-1', region=region, distributionType=HYDROSTATIC, 
    field='', magnitude=(h_up*5.2/1000./144.), hZero=botPoint+(h_up), hReference = botPoint)

x_down = int(h_down/12.) # Downstream hydraulic pressure
if (x_down< 10.):
    loadString = 'Water0'+str(x_down)
else:
    loadString = 'Water'+str(x_down)
    
mdb.models[namemodel].Pressure(name='Downstream', 
    createStepName='Step-1', region=region, distributionType=HYDROSTATIC, 
    field='', magnitude=(-h_down*5.2/1000./144.), hZero=botPoint+(h_down), hReference = botPoint)

a.regenerate()



mdb.Job(name=namejob, model=namemodel, description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=cpuCount, 
    numDomains=cpuCount, numGPUs=0)

mdb.jobs[namejob].writeInput(consistencyChecking=OFF)
mdb.save()
mdb.close()
