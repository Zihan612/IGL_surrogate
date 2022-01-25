# -*- coding: utf-8 -*-
"""
Created on Wed Oct 30 14:56:05 2019

@author: RDCHLTBF
"""

import numpy as np
#from fatigueAlgorithms import*
from odbAccess import *
from abaqusConstants import *
from odbMaterial import *
from odbSection import *
import csv
import sys
import time

#######################################################################
# Submit global Abaqus job, required every Monte Carlo step
#######################################################################


#######################################################################
# Record global and local Abaqus job information. Required once before Monte Carlo run
#######################################################################

np.set_printoptions(linewidth=np.nan)

workPath = 'C:\\Temp'
odbName = 'Submodel_global.odb'
LocalJobName = 'local.inp'
LocalInpName = 'Submodel_local.inp' 
Step = 'Step-1'
LocalStep = 'Step-1'

SectionOffset        = [                             0,                              0]  # 0 means middle, -1 means bottom surface, 1 means top surface
SectionThickness = [                          0.5,                           0.5]
SectionNormal      = [                    [0,0,1],                    [-1,0,0] ] 
partNames            = [             'PART-1-1',             'PART-1-1'  ]
setNames             = ['HORBOUNDARY', 'VERTBOUNDARY']  # The set names are assumed to be the same in the global and local .cae
setNames=[x.upper() for x in setNames] #Capitalize all GaugeNames, All sets in Abaqus are in ALL CAPS
localInstanceName =  'Part-1-1' # as it appears in .inp file
ActiveSteps = ['Step-1']

##################################################################
# Read coordinates, displacements, and rotations of all nodes on global boundary. Required every Monte Carlo step (not coordinates though)
##################################################################

odb = openOdb(workPath+'\\'+ odbName, readOnly=FALSE) # Open the Odb File
Frame=odb.steps[ActiveSteps[0]].frames[-1]  # There is only one step

print('\n\n\n\n')

# Get all nodes in setNames for global
GlobalNodes = {}
for i,iSet in enumerate(setNames):
    GlobalNodes[iSet] = {}
    GlobalNodes[iSet]['BoundaryNodes'] = odb.rootAssembly.nodeSets[iSet].nodes[0]  # List of all nodes along local boundary

# Get coordinates, displacements, and rotations
DisplField = Frame.fieldOutputs
RotField    = Frame.fieldOutputs['UR']
DisplField = Frame.fieldOutputs['U']

for key,values in GlobalNodes.items():
    iSet = values['BoundaryNodes']
    GlobalNodes[key]['Coords'] = np.zeros((len(iSet),3))
    GlobalNodes[key]['Displ']     = np.zeros((len(iSet),3))
    GlobalNodes[key]['Rot']       = np.zeros((len(iSet),3))
    for j,jNode in enumerate(iSet):
        GlobalNodes[key]['Coords'] [j,:]   = jNode.coordinates
        GlobalNodes[key]['Displ'] [j,:]       = DisplField.getSubset(region=jNode).values[0].data
        GlobalNodes[key]['Rot'] [j,:]         = RotField.getSubset(region=jNode).values[0].data

odb.close()

for key,values in GlobalNodes.items():
    del GlobalNodes[key]['BoundaryNodes'] 


##################################################################
# Find top and bottom surface displacements. Required every Monte Carlo step
##################################################################

# Provide vector rotation transformations
def xRot(CartesianVector,theta):
    Rx = np.array([ [ 1,                    0,                    0],
                              [0, np.cos(theta), -np.sin(theta)],
                              [0, np.sin(theta),  np.cos(theta)]])
    return np.dot(Rx, CartesianVector)

def yRot(CartesianVector,theta):
    Ry = np.array([ [ np.cos(theta), 0, np.sin(theta)],
                              [                    0, 1,                   0],
                              [-np.sin(theta), 0, np.cos(theta)]])
    return np.dot(Ry, CartesianVector)

def zRot(CartesianVector,theta):
    Rz = np.array([ [ np.cos(theta), -np.sin(theta), 0],
                              [  np.sin(theta), np.cos(theta), 0],
                              [                     0,                    0, 1]])
    return np.dot(Rz, CartesianVector)

DisplTop = []
CoordsTop = []
DisplBot = []
CoordsBot = []
for key,values in GlobalNodes.items():
    i = setNames.index(key)
    iSet = values['Displ']
    # Vector from shell node to top surface
    distTop = SectionThickness[i]*0.5*-(-1+SectionOffset[i])*np.linalg.norm(np.array(SectionNormal[i])) 
    dirTop = np.array(SectionNormal[i])
    # Vector from shell node to bottom surface
    distBot = SectionThickness[i]*0.5*(1+SectionOffset[i])*np.linalg.norm(np.array(SectionNormal[i])) 
    dirBot = -np.array(SectionNormal[i])
    
    # Update top and bottom surface coordinates for set
    Coords_temp = values['Coords']
    GlobalNodes[key]['CoordsTop'] = Coords_temp + distTop*dirTop
    GlobalNodes[key]['CoordsBot'] = Coords_temp + distBot*dirBot
    
    iRotVals = values['Rot']
    DisplTop_temp = np.zeros((np.shape(iSet)[0],3))
    DisplBot_temp = np.zeros((np.shape(iSet)[0],3))
    
    for j,jNode in enumerate(iSet):
        # Rotate vector to top surface
        vTop = distTop*dirTop
        vTop = xRot(vTop,iRotVals[j,0])
        vTop = yRot(vTop,iRotVals[j,1])
        vTop = zRot(vTop,iRotVals[j,2])
        
        DisplTop_temp[j,:] = jNode + (vTop - distTop*dirTop) # Add rotational top surface displacement to shell node displacement
        
        # Rotate vector to bottom surface
        vBot = distBot*dirBot
        vBot = xRot(vBot,iRotVals[j,0])
        vBot = yRot(vBot,iRotVals[j,1])
        vBot = zRot(vBot,iRotVals[j,2])
        
        DisplBot_temp[j,:] = jNode + (vBot - distBot*dirBot) # Add rotational bottom surface displacement to shell node displacement
    
    GlobalNodes[key]['DisplTop'] =  DisplTop_temp
    GlobalNodes[key]['DisplBot'] =   DisplBot_temp
    
print('\n Coords = \n\n\n')
print(GlobalNodes['HORBOUNDARY']['Coords'])
print('\n CoordsTop = \n\n\n')
print(GlobalNodes['HORBOUNDARY']['CoordsTop'])

##################################################################
# Read local .inp file to find all nodes in setNames. Required once before Monte Carlo run
##################################################################

def comma2intlist(string):
        tmp = [y.strip() for y in string.split(',')]
        tmp = filter(None,tmp)
        return [int(y) for y in tmp]

def readSet(LineString,generateFlag):
        Set = []
        if generateFlag == True: #The nodes are not explicitly listed but the pattern is given
            tmp = comma2intlist(x)
            Set = list(np.arange(tmp[0],tmp[1]+1,tmp[2]))
        else:
            Set += comma2intlist(x)
        return Set

setFlags = [False for i in setNames]
generateFlag = False

LocalNodes = {}
for x in setNames:
    LocalNodes[x] = {}
    LocalNodes[x]['Index'] = []
    LocalNodes[x]['Coords'] = []
LocalNodesFlat = []

with open(workPath+'\\'+LocalInpName,'r') as f: 
    for index, x in enumerate(f):
        if '*' in x:
            setFlags = [False for i in setNames]
            generateFlag = False
        
        for indy,y in enumerate(setNames):
            if setFlags[indy] == True:
                nodes = readSet(x,generateFlag)
                LocalNodes[y]['Index'] += nodes
                LocalNodesFlat += nodes
            
            if '*Nset, nset='+y in x:
                setFlags[indy] = True
                if 'generate' in x:
                    generateFlag = True
f.close()
LocalNodesFlat = sorted(set(LocalNodesFlat))

###############################################################
# Get coordinates of all nodes in local domain. Required once before Monte Carlo run
###############################################################

def comma2floatlist(string):
        tmp = [y.strip() for y in string.split(',')]
        tmp = filter(None,tmp)
        return [int(tmp[0]), float(tmp[1]), float(tmp[2]), float(tmp[3])]

# Assumed there is only one instance, xyz coordinates,
LocalNodesFlatCoords = np.empty((len(LocalNodesFlat),3))

LocalNodesFlatIndex = 0
NodeCoords = False
with open(workPath+'\\'+LocalInpName,'r') as f: 
    for index, x in enumerate(f):
     
        if '*' in x: 
            NodeCoords = False
        if '*Node' in x:
            NodeCoords=True
            continue # No coordinate information on this line
            
        if NodeCoords == True:
            NodeInfo = comma2floatlist(x)
            if NodeInfo[0] == LocalNodesFlat[LocalNodesFlatIndex]:
                LocalNodesFlatCoords[LocalNodesFlatIndex] = NodeInfo[1:]
                LocalNodesFlatIndex = LocalNodesFlatIndex+1
                if LocalNodesFlatIndex > len(LocalNodesFlat)-1:
                    break
f.close()
print(LocalNodesFlatCoords)

# Adjust coordinates by translation and rotation in assembly
flagTranslationInfo = False
translations = np.zeros(3)
with open(workPath+'\\'+LocalInpName,'r') as f: 
    for index, x in enumerate(f):
        if '*Instance, name='+localInstanceName in x:
            flagTranslationInfo = True
            continue
        if flagTranslationInfo == True:
            tmp = [y.strip() for y in x.split(',')]
            tmp = filter(None,tmp)
            translations = np.array([float(tmp[0]), float(tmp[1]), float(tmp[2])])
            break
f.close()

for i,iCoord in enumerate(LocalNodesFlatCoords):
    LocalNodesFlatCoords[i,:] = iCoord + translations

# Match coordinates to sets
for key, value in LocalNodes.items():
    LocalNodes[key]['Coords'] = np.zeros((len(value['Index']),3))
    FlatInd = 0
    for i,iNode in enumerate(value['Index']):
        for j,jNode in enumerate(LocalNodesFlat[FlatInd:]):
            if iNode == jNode:
                LocalNodes[key]['Coords'][i,:] = LocalNodesFlatCoords[j+FlatInd,:]
                FlatInd=j


#################################################
# Interpolate between displacement values in global domain. Required every Monte Carlo step
#################################################

from scipy.interpolate import griddata
from scipy.interpolate import Rbf

for key, value in GlobalNodes.items():
    points = np.concatenate( (value['CoordsTop'], value['Coords'], value['CoordsBot'] ), axis = 0)
    values = np.concatenate( (value['DisplTop'], value['Displ'], value['DisplBot'] ), axis = 0)
    
    grid_x = LocalNodes[key]['Coords'][:,0]
    grid_y = LocalNodes[key]['Coords'][:,1]
    grid_z = LocalNodes[key]['Coords'][:,2] 
    
    # RBF used to interpolate. Due to computer precision issues, nodes in the local domain may not fall within the 
    # convex hull of the nodes in the global domain. Thus strict interpolation can't be used. Consequently, reduced
    # basis function (rbf) splines are used that interpolate well but can also provided limited extrapolation as we need.
    
    rbfx = Rbf(points[:,0], points[:,1], points[:,2], values[:,0])
    dx = rbfx(grid_x, grid_y, grid_z)
    rbfy = Rbf(points[:,0], points[:,1], points[:,2], values[:,1])
    dy = rbfy(grid_x, grid_y, grid_z)
    rbfz = Rbf(points[:,0], points[:,1], points[:,2], values[:,2])
    dz = rbfz(grid_x, grid_y, grid_z)
    
    LocalNodes[key]['Displ'] = np.stack([dx,dy,dz],axis=1)
    print(LocalNodes[key]['Displ'])
    print(dx[-10:])
    print(dy[-10:])
    print(dz[-10:])
    print(LocalNodes[key]['Coords'][-10:,:])

    
##############################################
# Apply displacement values to .inp file. Required every Monte Carlo step
##############################################

def getlines(fobj,line1):
     for line in iter(fobj.readline,''):  #This is necessary to get `fobj.tell` to work
         yield line
         if line == line1:
            #pos = fobj.tell()
            #next_line = next(fobj)
            #fobj.seek(pos)
            print('Found it')
            break

with open(workPath+'\\'+LocalInpName) as fin, open(LocalJobName,'w') as fout:     
    BCCounter = 5
    for line in fin:
            fout.write(line)
            if line == '*End Instance\n':  # Create all of the nodesets for the boundary conditions
                for key, value in LocalNodes.items():
                    for NodeNum in value['Index']: 
                        fout.write('*Nset, nset=GL-'+str(NodeNum)+', instance='+localInstanceName+'\n')
                        fout.write(str(NodeNum)+',\n')
            
            if line == '*Step, name='+ LocalStep +', nlgeom=NO\n': # Create boundary conditions after step definition
                BCCounter = 0
            BCCounter = BCCounter + 1
            if BCCounter == 4:
                fout.write('** BOUNDARY CONDITIONS\n')
                for key, value in LocalNodes.items():
                    for n, NodeNum in enumerate(value['Index']): 
                        fout.write('*Boundary\n')
                        fout.write('GL-'+str(NodeNum)+', 1, 1, '+ str(value['Displ'][n,0]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 2, 2, '+ str(value['Displ'][n,1]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 3, 3, '+ str(value['Displ'][n,2]) +'\n')
f.close()

#############################################
# Run local .inp file. Required every Monte Carlo step
#############################################