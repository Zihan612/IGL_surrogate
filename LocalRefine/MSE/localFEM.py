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


index = int(sys.argv[-4])
# index=1

cpuCount = int(sys.argv[-3])

workPath = sys.argv[-2]
# workPath = 'C:\ZIHAN\GreenupCrackEstimation\A1201\LocalRefine\MSE\Files\Model1'
# workPath = workPath + "\\S" +str(index)+ "\\"
os.chdir(workPath)

crack_length = float(sys.argv[-1])

namejob = 'GreenupLocal'+str(index)
mdb=openMdb(namejob+ '.cae')

namemodel = 'LocModel-'+str(index)
mdb.models.changeKey(fromName='Model-1', toName=namemodel)

p = mdb.models[namemodel].parts['Crack']
s = p.features['Shell planar-1'].sketch
mdb.models[namemodel].ConstrainedSketch(name='__edit__', objectToCopy=s)
s1 = mdb.models[namemodel].sketches['__edit__']
g, v, d, c = s1.geometry, s1.vertices, s1.dimensions, s1.constraints
s1.setPrimaryObject(option=SUPERIMPOSE)
p.projectReferencesOntoSketch(sketch=s1, 
    upToFeature=p.features['Shell planar-1'], filter=COPLANAR_EDGES)

s1.FixedConstraint(entity=g.findAt((23.75, 10.875)))

d[1].setValues(value=crack_length, )
s1.unsetPrimaryObject()

p.features['Shell planar-1'].setValues(sketch=s1)
del mdb.models[namemodel].sketches['__edit__']

p.regenerate()
a = mdb.models[namemodel].rootAssembly
a.regenerate()

mdb.save()

mdb.Job(name=namejob, model=namemodel, description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=cpuCount, 
    numDomains=cpuCount, numGPUs=0)

mdb.jobs[namejob].writeInput(consistencyChecking=OFF)
mdb.close()