# -*- coding: utf-8 -*-
"""
Created on Fri Dec 18 09:56:49 2020

@author: RDCHLTBF
"""

import subprocess
import os
import sys
import math
import numpy as np
import shutil
        
def createSubStructureInput(workPath,scriptPath,abaqusVersion,caeName,modelName,stepName,instanceName, partNames,
                            boundary,cpuCount,SSSavedFlag):
    print('CreateSubStructureInput !!!!!!!!####################@@@@@@@@@@@@@@@@@')
    print(caeName)
    open(workPath+'\\stderr.txt', 'w').close() #erase contents of error file
    with open(workPath+'\\stdout.txt','wb') as out, open(workPath+'\\stderr.txt','wb') as err:
        abaqusCall1 = 'cmd.exe /c '+abaqusVersion+' cae noGUI=PrepareCAEStaticCond.py' \
                        +' -- '+abaqusVersion+' '+workPath+' '+scriptPath+' ' \
                        +modelName+' '+stepName+' '+boundary+' ' \
                        +str(cpuCount)+' '+caeName
        process1 = subprocess.Popen(abaqusCall1, cwd=scriptPath,stdout=out,stderr=err)
        process1.wait()
        process1.terminate()
        try:
            process1.kill()
        except:
            pass
        del process1
    if os.stat(workPath+'\\stdout.txt').st_size != 0:
        with open(workPath+'\\'+'stdout.txt') as f:
            for cnt, line in enumerate(f):
                print(line)
        with open(workPath+'\\'+'stderr.txt') as f:
            for cnt, line in enumerate(f):
                print(line)
        sys.exit('Abaqus Error! See above')
        
    jobName = caeName[:-4]+'_mat'+'.inp'
    os.chdir(workPath)   
    if SSSavedFlag == False:
        # Modify .inp file: add line to generate substructured matrix
        lookStars = False
        with open(workPath+caeName[:-4]+'_Sub.inp') as fin, open(workPath+caeName[:-4]+'_mat.inp','w') as fout:     
            for line in fin:
                fout.write(line)
                if '*Step, name='+stepName in line: 
                    lookStars = True
                if lookStars == True:
                    if 'Substructure' in line:
                        lookStars = False
                        fout.write('*SUBSTRUCTURE MATRIX OUTPUT, OUTPUT FILE=USER DEFINED, FILE NAME='+caeName[:-4]+'_mat, SLOAD=YES,STIFFNESS=YES\n')
        
        # Submit job
        ############         
    
        open(workPath+'\stderr.txt', 'w').close() #erase contents of error file
        with open(workPath+'\stdout.txt','wb') as out, open(workPath+'\stderr.txt','wb') as err:
            abaqusCall1 = 'cmd.exe /c '+abaqusVersion+' cae noGUI='+scriptPath+'\\SubmitAbaqusJob.py' \
                            +' -- '+workPath+' '+jobName
            process1 = subprocess.Popen(abaqusCall1, cwd=workPath,stdout=out,stderr=err)
            process1.wait()
            process1.terminate()
            try:
                process1.kill()
            except:
                pass
            del process1
            
        if os.stat(workPath+'\stdout.txt').st_size != 0:
            with open(workPath+'\\'+'stdout.txt') as f:
                for cnt, line in enumerate(f):
                    print(line)
            with open(workPath+'\\'+'stderr.txt') as f:
                for cnt, line in enumerate(f):
                    print(line)
            sys.exit('Abaqus Error! See above')
        
    else:
        pass # The .mtx file already exists and contains the necessary information
    ############################################################################
    #       Get stiffness matrix, load matrix, and nodal mapping from uncondensed to          *********
    #        condensed nodes from the .mtx file                                                                       *********
    ############################################################################
    
    FNAME = jobName[:-4]+'.mtx'
    
    atStiff = False
    atLoad = False
    recordNodesDOFs = True
    
    nNodeLines = 0
    StiffMat1D = []
    LoadName = ''
    LoadMat = []
    LoadMatNames = ['Load_0']
    condNodesIndex = []
    DOFs = []
    nodesInd_inp = []
    with open(workPath+FNAME,'r') as f: 
        for index, x in enumerate(f):
            x = x.rstrip()
            # the third row has the number of nodes
            if index == 2:
                nNodes = int(x.split(',')[1].split()[-1])
                nNodeLines = math.ceil(nNodes/10)
                
            # there are 10 nodes on each line starting on line 5
            elif 4 <= index and index < 4+nNodeLines:
                tmp = [y.strip() for y in x[2:].split(',')]
                tmp = filter(None, tmp)
                nodesInd_inp += [int(y) for y in tmp] # This is a list of the index of all uncondensed nodes
            # Read stiffness matrix info in
            elif index >= 4+nNodeLines:
                if '*MATRIX' in x:
                    atStiff = True
                    continue
                elif '** SUBSTRUCTURE' in x and atLoad == False:
                    atStiff = False
                    atLoad = True
                    
                    LoadName = x.split()[-1].strip()
                    LoadMatNames += [LoadName]
                    LoadMat += [np.zeros(6*len(nodesInd_inp))]
                    continue
                elif '** SUBSTRUCTURE' in x and atLoad == True:
                    recordNodesDOFs = False
                    LoadName = x.split()[-1].strip()
                    LoadMatNames += [LoadName]
                    LoadMat += [np.zeros(6*len(nodesInd_inp))]
                    continue
                if '***CLOAD' in x:
                    continue
                
                if atStiff == True:
                    tmp = [y.strip() for y in x.split(',')]
                    tmp = filter(None,tmp)
                    StiffMat1D += [float(y) for y in tmp]
                    
                elif atLoad == True:
                    
                    tmp = [y.strip() for y in x[2:].split(',')]
                    LoadMat[-1][(int(tmp[0])-1)*6+int(tmp[1])-1] = float(tmp[2])
                    
                    
                    if recordNodesDOFs == True:
                        condNodesIndex += [int(tmp[0])]
                        DOFs += [int(tmp[1])]
    LoadMat.insert(0, np.zeros(6*len(nodesInd_inp)))
    # Add dof information to stiffness matrix
    # Create full stiffness matrix from 1d vector
    tmp = []
    for i, iLoad in enumerate(condNodesIndex):
        tmp += [nodesInd_inp[iLoad-1]]
    condNodesIndex = tmp
    
    Stiff_len = int((-1 + np.sqrt(1 + 8*len(StiffMat1D)))/2) #Inverse of 1+2+...+n = n*(n+1)/2
    StiffMat = np.zeros((Stiff_len,Stiff_len))
    colEnd = 1
    col = 0
    row = 0
    for i, val in enumerate(StiffMat1D):
        StiffMat[row,col] = val
        StiffMat[col,row] = val
        col += 1
        if col == colEnd:
            row += 1
            colEnd += 1
            col = 0    
    
    ##################################################################
    # Get condensed coordinates for nodesInd_inp
    ##################################################################
    
    def comma2floatlist(string):
                tmp = [y.strip() for y in string.split(',')]
                #tmp = filter(None,tmp)
                return [int(tmp[0]), float(tmp[1]), float(tmp[2]), float(tmp[3])]
        
    # Assumed there is only one instance, xyz coordinates,
    nodesCoord_inp = np.empty((len(nodesInd_inp),3))
    
    nodesIndex = 0
    nodeCoordsFlag = False
    partFlag = False
    transAxisFlag = False
    AxisFlag = False
    translations = np.zeros(3)
    x1Axis = np.array([1.,0.,0.])
    y1Axis = np.array([0.,1.,0.])
    inpName = jobName
    with open(workPath+'\\'+inpName,'r') as f: 
        for index, x in enumerate(f):
            
            if '** PART INSTANCE: '+instanceName in x:
                partFlag = True
            if '*System' in x and partFlag == True:
                transAxisFlag = True
                continue
            elif '*' in x: 
                nodeCoordsFlag = False
            if '*Node' in x:
                AxisFlag = False
                nodeCoordsFlag = True
                continue # No coordinate information on this line
                
            if nodeCoordsFlag == True and partFlag == True and AxisFlag == False and transAxisFlag == False:
                nodeInfo = comma2floatlist(x)
                if nodeInfo[0] == nodesInd_inp[nodesIndex]:
                    nodesCoord_inp[nodesIndex] = nodeInfo[1:]
                    nodesIndex = nodesIndex+1
                    if nodesIndex == len(nodesInd_inp):
                        break
            if transAxisFlag == True:
                tmp = [y.strip() for y in x.split(',')]
                translations = np.array([float(tmp[0]), float(tmp[1]), float(tmp[2])])
                x1Axis = np.array([float(tmp[3]), float(tmp[4]), float(tmp[5])])
                transAxisFlag = False
                y1Axis = y1Axis + translations
                AxisFlag = True
                continue
            
            if AxisFlag == True:
                tmp = [y.strip() for y in x.split(',')]
                y1Axis = np.array([float(tmp[0]), float(tmp[1]), float(tmp[2])])
                AxisFlag = False
                
    f.close()
    
    x = np.array([1.,0.,0.])
    y = np.array([0.,1.,0.])
    u = x1Axis - translations
    v = y1Axis - translations
    
    angleNormal = np.arccos(np.dot(x,u))
    crossNormal = np.cross(x,u)
    if np.linalg.norm(crossNormal) > 10e-8:
        crossNormal = crossNormal/np.linalg.norm(crossNormal)
    
    I = np.identity(3)
    CNormal = np.array([[0,-crossNormal[2],crossNormal[1]],
                  [crossNormal[2],0,-crossNormal[0]],
                  [-crossNormal[1],crossNormal[0],0]])
    RNormal = I + CNormal*np.sin(angleNormal) + np.dot(CNormal,CNormal)*(1-np.cos(angleNormal))
    y = np.dot(RNormal,y)
    angleUp = np.arccos(np.dot(y,v))
    crossUp = np.cross(y,v)
    if np.linalg.norm(crossUp) > 10e-8:
        crossUp = crossUp/np.linalg.norm(crossUp)

    CUp = np.array([[0,-crossUp[2],crossUp[1]],
                  [crossUp[2],0,-crossUp[0]],
                  [-crossUp[1],crossUp[0],0]])
    RUp = I + CUp*np.sin(angleUp) + np.dot(CUp,CUp)*(1-np.cos(angleUp))
    
    R = np.dot(RUp,RNormal)
    
    for i,iCoord in enumerate(nodesCoord_inp):
        nodesCoord_inp[i,:] = np.dot(R,iCoord.T)
    
    for i,iCoord in enumerate(nodesCoord_inp):
        nodesCoord_inp[i,:] = iCoord + translations
        
    #os.remove(os.path.join(self.workPath, FNAME))
    np.save(workPath+jobName[:-4]+'LoadMat',np.array(LoadMat))
    np.save(workPath+jobName[:-4]+'LoadMatNames',np.array(LoadMatNames))
    np.save(workPath+jobName[:-4]+'StiffMat',StiffMat)
    np.save(workPath+jobName[:-4]+'nodesInd_inp',np.array(nodesInd_inp))
    np.save(workPath+jobName[:-4]+'nodesCoord_inp',np.array(nodesCoord_inp))
    np.save(workPath+jobName[:-4]+'condNodesIndex',np.array(condNodesIndex))
    np.save(workPath+jobName[:-4]+'DOFs',np.array(DOFs))
    ##################################################################
    # Read uncondensed .inp file to find all nodes in boundary. Required once before Monte Carlo run
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
    
    setFlags = False
    generateFlag = False
    
    nodesInd_inp2 = []
    NSCInpName = caeName[:-4]+'.inp'
    with open(workPath+'\\'+NSCInpName,'r') as f: 
        for index, x in enumerate(f):
            if '*' in x:
                setFlags = False
                generateFlag = False
            
            if setFlags == True:
                nodes = readSet(x,generateFlag)
                nodesInd_inp2 += nodes
            
            if '*Nset, nset='+boundary.upper() in x:
                setFlags = True
                if 'generate' in x:
                    generateFlag = True
    f.close()
    nodesInd_inp2 = sorted(set(nodesInd_inp2))
    if len(nodesInd_inp2) == 0:
        sys.exit('No nodes found for input boundary: '+boundary)
    np.save(workPath+caeName[:-4]+'nodesInd_inp2',np.array(nodesInd_inp2))
    
    ###############################################################
    # Get coordinates of all uncodensed nodes
    ###############################################################
    
    def comma2floatlist(string):
            tmp = [y.strip() for y in string.split(',')]
            #tmp = filter(None,tmp)
            return [int(tmp[0]), float(tmp[1]), float(tmp[2]), float(tmp[3])]
    
    # Assumed there is only one instance, xyz coordinates,
    nodesCoord_inp2 = np.empty((len(nodesInd_inp2),3))
    
    nodesIndex = 0
    nodeCoordsFlag = False
    partFlag = False
    with open(workPath+'\\'+NSCInpName,'r') as f: 
        for index, x in enumerate(f):
         
            if '*' in x: 
                nodeCoordsFlag = False
            if '*Part, name='+partNames in x:
                partFlag = True
            if '*Node' in x:
                nodeCoordsFlag = True
                continue # No coordinate information on this line
                
            if nodeCoordsFlag == True and partFlag == True:
                nodeInfo = comma2floatlist(x)
                if nodeInfo[0] == nodesInd_inp2[nodesIndex]:
                    nodesCoord_inp2[nodesIndex] = nodeInfo[1:]
                    nodesIndex = nodesIndex+1
                    if nodesIndex == len(nodesInd_inp2):
                        break
    f.close()
    
    # Adjust coordinates by translation and rotation in assembly
    flagTranslationInfo = False
    flagRotationInfo = False
    instanceFound = False
    translations = np.zeros(3)
    rotations = np.zeros(7)
    with open(workPath+'\\'+NSCInpName,'r') as f: 
        for index, x in enumerate(f):
            if '*' in x:
                flagTranslationInfo = False
                flagRotationInfo = False
                if '*Instance, name='+instanceName in x:
                    flagTranslationInfo = True
                    instanceFound = True
                    continue
            if flagTranslationInfo == True and '*' not in x:
                tmp = [y.strip() for y in x.split(',')]
                #tmp = filter(None,tmp)
                translations = np.array([float(tmp[0]), float(tmp[1]), float(tmp[2])])
                flagTranslationInfo = False
                flagRotationInfo = True
                continue
            if flagRotationInfo == True and '*' not in x:
                tmp = [y.strip() for y in x.split(',')]
                #tmp = filter(None,tmp)
                rotations = np.array([float(tmp[0]), float(tmp[1]), float(tmp[2]), float(tmp[3]), float(tmp[4]), float(tmp[5]), float(tmp[6])])
                flagRotationInfo = False
                break
        if instanceFound == False:
            sys.exit('globinstance '+instanceName+' not found in '+workPath+'\\'+NSCInpName)
                
    f.close()

    for i,iCoord in enumerate(nodesCoord_inp2):
        nodesCoord_inp2[i,:] = iCoord + translations

    v = rotations[3:6]-rotations[0:3]
    theta = rotations[6]
    I = np.identity(3)
    C = np.array([[0,-v[2],v[1]],
                  [v[2],0,-v[0]],
                  [-v[1],v[0],0]])
    R = I + C*np.sin(np.radians(theta)) + np.dot(C,C)*(1-np.cos(np.radians(theta)))
    
    for i,iCoord in enumerate(nodesCoord_inp2):
        nodesCoord_inp2[i,:] = np.dot(R,iCoord.T).T 
    np.save(workPath+caeName[:-4]+'nodesCoord_inp2',nodesCoord_inp2)
    
def modify_local_crack(info):
    shutil.copy(info['scriptPath']+info['locInfo']['baseCAEName'], info['locPath']+info['locInfo']['CAEName'])
    print('Modify local crack')
    normalStr = ''.join(np.array2string(info['locInfo']['normal'], separator=',').split())
    propStr = ''.join(np.array2string(info['locInfo']['propDir'], separator=',').split())
    translStr = ''.join(np.array2string(info['locInfo']['translation'], separator=',').split())
    open(info['workPath']+'\\stderr.txt', 'w').close() #erase contents of error file
    with open(info['locPath']+'stdout.txt','wb') as out, open(info['locPath']+'stderr.txt','wb') as err:
        abaqusCall1 = 'cmd.exe /c '+info['IGLInfo']['abaqusVersion']+' cae noGUI=ModifyLocalCrack.py' \
                        +' -- '+info['IGLInfo']['abaqusVersion']+' '+info['locPath']+' '\
                        +info['locInfo']['baseCAEName']+' '+str(info['identifier'])+' ' \
                        +str(info['locInfo']['crackLength'])+' '+str(info['locInfo']['thickness'])+' ' \
                        +normalStr+' '+propStr+' '+translStr
        process1 = subprocess.Popen(abaqusCall1, cwd=info['scriptPath'],stdout=out,stderr=err)
        process1.wait()
        process1.terminate()
        try:
            process1.kill()
        except:
            pass
        del process1
    if os.stat(info['locPath']+'\\stdout.txt').st_size != 0:
        with open(info['locPath']+'\\'+'stdout.txt') as f:
            for cnt, line in enumerate(f):
                print(line)
        with open(info['locPath']+'\\'+'stderr.txt') as f:
            for cnt, line in enumerate(f):
                print(line)
        sys.exit('Abaqus Error! See above')
    
    #os.remove(info['scriptPath']+info['locInfo']['CAEName']) 