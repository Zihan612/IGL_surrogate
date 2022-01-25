# -*- coding: utf-8 -*-
"""
Created on Wed Sep  2 14:45:55 2020

@author: RDCHLTBF
"""
import subprocess
import os
import time
import sys
import numpy as np

class ShellAbaqus:
    def __init__(self,abaqusVersion,workPath,scriptPath,info):
        self.workPath  = workPath
        self.scriptPath  = scriptPath
        self.saveResults = info['saveResults']
        self.instanceNames = info['instanceName']
        self.partNames = info['partNames']
        self.boundary  = info['boundary']
        self.globLocDomain = info['globLocDomain']
        self.stepName  = info['stepName']
        self.inpName = info['inpName']
        self.abaqusVersion = abaqusVersion
        
        self.odbName   = self.inpName[:-4]+'.odb'
        self.nodesInd_inp = []
        self.nodesCoord_inp = []
        self.nodesCoord_odb = []
        
        self.fileResultNames = info['inpName'][:-4]+'Result.txt'
        self.fileOutputName = info['inpName'][:-4]+'Output'
        
        self.get_nodes_input()
        self.get_nodes_coord_input()
        
    def printShell(self,text):
        with open(self.workPath+'stdout.txt', 'a') as f:
            print(text, file=f)
    
    def get_nodes_input(self):
        ##################################################################
        # Read local .inp file to find all nodes in boundary. Required once before Monte Carlo run
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
        
        self.nodesInd_inp = []
        self.nodesCoord_inp = []
        
        with open(self.workPath+'\\'+self.inpName,'r') as f: 
            for index, x in enumerate(f):
                if '*' in x:
                    setFlags = False
                    generateFlag = False
                
                if setFlags == True:
                    nodes = readSet(x,generateFlag)
                    self.nodesInd_inp += nodes
                
                if '*Nset, nset='+self.boundary.upper() in x:
                    setFlags = True
                    if 'generate' in x:
                        generateFlag = True
        f.close()
        self.nodesInd_inp = sorted(set(self.nodesInd_inp))
        self.printShell('Done finding local node indices')
        if len(self.nodesInd_inp) == 0:
            self.printShell('No nodes found for input boundary: '+self.boundary)
            sys.exit()
        
    
    def get_nodes_coord_input(self):
        ###############################################################
        # Get coordinates of all nodes in local domain. Required once before Monte Carlo run
        ###############################################################
        
        def comma2floatlist(string):
                tmp = [y.strip() for y in string.split(',')]
                #tmp = filter(None,tmp)
                return [int(tmp[0]), float(tmp[1]), float(tmp[2]), float(tmp[3])]
        
        # Assumed there is only one instance, xyz coordinates,
        self.nodesCoord_inp = np.empty((len(self.nodesInd_inp),3))
        
        nodesIndex = 0
        nodeCoordsFlag = False
        partFlag = False
        with open(self.workPath+'\\'+self.inpName,'r') as f: 
            for index, x in enumerate(f):
             
                if '*' in x: 
                    nodeCoordsFlag = False
                if '*Part, name='+self.partNames in x:
                    partFlag = True
                if '*Node' in x:
                    nodeCoordsFlag = True
                    continue # No coordinate information on this line
                    
                if nodeCoordsFlag == True and partFlag == True:
                    nodeInfo = comma2floatlist(x)
                    if nodeInfo[0] == self.nodesInd_inp[nodesIndex]:
                        self.nodesCoord_inp[nodesIndex] = nodeInfo[1:]
                        nodesIndex = nodesIndex+1
                        if nodesIndex > len(self.nodesInd_inp)-1:
                            break
        f.close()
        
        # Adjust coordinates by translation and rotation in assembly
        flagTranslationInfo = False
        flagRotationInfo = False
        instanceFound = False
        translations = np.zeros(3)
        rotations = np.zeros(7)
        with open(self.workPath+'\\'+self.inpName,'r') as f: 
            for index, x in enumerate(f):
                if '*' in x:
                    flagTranslationInfo = False
                    flagRotationInfo = False
                    if '*Instance, name='+self.instanceNames in x:
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
                self.printShell('globinstance '+self.instanceNames+' not found in '+self.workPath+'\\'+self.inpName)
                sys.exit()
                    
        f.close()

        for i,iCoord in enumerate(self.nodesCoord_inp):
            self.nodesCoord_inp[i,:] = iCoord + translations

        v = rotations[3:6]-rotations[0:3]
        theta = rotations[6]
        I = np.identity(3)
        C = np.array([[0,-v[2],v[1]],
                      [v[2],0,-v[0]],
                      [-v[1],v[0],0]])
        R = I + C*np.sin(np.radians(theta)) + np.dot(C,C)*(1-np.cos(np.radians(theta)))
        
        for i,iCoord in enumerate(self.nodesCoord_inp):
            self.nodesCoord_inp[i,:] = np.dot(R,iCoord.T).T
            
        np.save(self.workPath+'Coord_inp.npy',self.nodesCoord_inp)
        
        self.printShell('Done finding local node coordinates')
            
    def apply_displacement(self, AppliedDispl,iteration):        
        
    ##############################################
    # Apply displacement values to .inp file
    ##############################################  
        jobName = self.odbName[:-4]+'_job'+str(iteration)+'.inp'
        fieldFlag = False
        with open(self.workPath+'\\'+self.inpName) as fin, open(self.workPath+'\\'+jobName,'w') as fout:     
            BCCounter = 5
            for line in fin:
                if '**FIELD OUTPUT' in line:
                    fieldFlag = True
                if fieldFlag == True:
                    if '*End Step' in line:
                        fout.write('*Output, field\n')
                        fout.write('*Element Output, elset='+self.boundary+', directions=YES\n')
                        fout.write('NFORC\n')
                        fout.write('*End Step\n')
                    continue
                fout.write(line)
                if line == '*End Instance\n':  # Create all of the nodesets for the boundary conditions
                    for NodeNum in self.nodesInd_inp:
                        fout.write('*Nset, nset=GL-'+str(NodeNum)+', instance='+self.instanceNames+'\n')
                        fout.write(str(NodeNum)+',\n')
                
                if line == '*Step, name='+ self.stepName +', nlgeom=NO\n': # Create boundary conditions after step definition
                    BCCounter = 0
                BCCounter = BCCounter + 1
                if BCCounter == 4:
                    fout.write('** BOUNDARY CONDITIONS\n')
                    for n,NodeNum in enumerate(self.nodesInd_inp):
                        fout.write('*Boundary\n')
                        fout.write('GL-'+str(NodeNum)+', 1, 1, '+ str(AppliedDispl[n,0]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 2, 2, '+ str(AppliedDispl[n,1]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 3, 3, '+ str(AppliedDispl[n,2]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 4, 4, '+ str(AppliedDispl[n,3]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 5, 5, '+ str(AppliedDispl[n,4]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 6, 6, '+ str(AppliedDispl[n,5]) +'\n')
        fout.close()
            
    def apply_extra_load(self,ExtraLoad,iteration):
        
    ##############################################
    # Apply global extra load values to .inp file
    ############################################## 
        jobName = self.odbName[:-4]+'_job'+str(iteration)+'.inp'
        fieldFlag = False
        with open(self.workPath+'\\'+self.inpName) as fin, open(self.workPath+'\\'+jobName,'w') as fout:     
            for line in fin:
                if '**FIELD OUTPUT' in line:
                    fieldFlag = True
                if fieldFlag == True:
                    if '*End Step' in line:
                        fout.write('*Output, field\n')
                        fout.write('*Node Output, nset='+self.boundary+'\n')
                        fout.write('RF, RM, U\n')
                        fout.write('*End Step\n')
                    continue
                fout.write(line)
                if line == '*End Instance\n':  # Create all of the nodesets for the boundary conditions
                    for NodeNum in self.nodesInd_inp:
                        fout.write('*Nset, nset=GL-'+str(NodeNum)+', instance='+self.instanceNames+'\n')
                        fout.write(str(NodeNum)+',\n')
                
                if line == '** LOADS\n':
                    for n, NodeNum in enumerate(self.nodesInd_inp):
                        fout.write('*Cload\n')
                        fout.write('GL-'+str(NodeNum)+', 1, '+ str(ExtraLoad[n,0]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 2, '+ str(ExtraLoad[n,1]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 3, '+ str(ExtraLoad[n,2]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 4, '+ str(ExtraLoad[n,3]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 5, '+ str(ExtraLoad[n,4]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 6, '+ str(ExtraLoad[n,5]) +'\n')
        fout.close()
        
    def run_job(self,iteration):
        if iteration > 0 and self.saveResults == False:
            test = os.listdir(self.workPath)
            for item in test: 
                if item.endswith('job'+str(iteration-1)+'.odb'):
                    #os.remove(os.path.join(self.workPath, item))
                    pass
            
        jobName = self.odbName[:-4]+'_job'+str(iteration)
        os.chdir(self.workPath)
        self.printShell('Running job')
        try:
            subprocess.check_call(self.abaqusVersion+" job="+jobName+" cpus=2 output_precision=full", stdin=None, stdout=None, stderr=None, shell=True, timeout = 20000.0)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError) as e:
            pass
        
        time.sleep(5)
        while not os.path.exists(self.workPath+'\\'+jobName+'.sta'):
            time.sleep(1)
            if not os.path.exists(self.workPath+'\\'+jobName+'.lck'):
                self.printShell('Abaqus simulation exited without convergence')
                sys.exit()
    
        completedFlag = False
        while completedFlag == False:
            if os.path.isfile(self.workPath+'\\'+jobName+'.sta'):
                # read file
                with open(self.workPath+'\\'+jobName+'.sta','r') as f: 
                    self.printShell('checking .sta')
                    for index, x in enumerate(f):
                        if 'COMPLETED SUCCESSFULLY' in x:
                            completedFlag = True
                        elif 'NOT BEEN COMPLETED' in x:
                            self.printShell(jobName+' was aborted in Abaqus')
                            sys.exit()
                            
                f.close()
                time.sleep(1)
            else:
                raise ValueError("%s isn't a file!" % self.workPath+'\\'+jobName+'.sta')
        '''
        open(self.workPath+'\stderr.txt', 'w').close() #erase contents of error file
        with open(self.workPath+'\stdout.txt','wb') as out, open(self.workPath+'\stderr.txt','wb') as err:
            abaqusCall1 = 'cmd.exe /c '+self.abaqusVersion+' cae noGUI='+self.scriptPath+'\\SubmitAbaqusJob.py' \
                            +' -- '+self.workPath+' '+jobName
            process1 = subprocess.Popen(abaqusCall1, cwd=self.workPath,stdout=out,stderr=err)
            process1.wait()
            process1.terminate()
            del process1
            
        if os.stat(self.workPath+'\stdout.txt').st_size != 0:
            with open(self.workPath+'\\'+'stdout.txt') as f:
                for cnt, line in enumerate(f):
                    print(line)
            with open(self.workPath+'\\'+'stderr.txt') as f:
                for cnt, line in enumerate(f):
                    print(line)
            sys.exit('Abaqus Error! See above')
        '''
            
    def GetOdbResults(self, iteration):
        t_odbName = self.odbName[:-4]+'_job'+str(iteration)+'.odb'
        open(self.workPath+'\\stderr.txt', 'w').close() #erase contents of error file
        with open(self.workPath+'\\stdFlag.txt','wb') as out, open(self.workPath+'\\stderr.txt','ab') as err:
            abaqusCall1 = 'cmd.exe /c '+self.abaqusVersion+' cae noGUI=GetOdbResults.py' \
                            +' -- '+self.fileResultNames+' '+self.fileOutputName+' '+self.boundary.upper()+' '+self.stepName+' ' \
                            +t_odbName+' '+self.workPath
            process1 = subprocess.Popen(abaqusCall1, cwd=self.scriptPath,stdout=out,stderr=err)
            process1.wait()
            process1.terminate()
            del process1
        if os.stat(self.workPath+'\\stdFlag.txt').st_size != 0:
            with open(self.workPath+'\\'+'stdFlag.txt') as f:
                for cnt, line in enumerate(f):
                    print(line)
            with open(self.workPath+'\\'+'stderr.txt') as f:
                for cnt, line in enumerate(f):
                    print(line)
            self.printShell('Abaqus Error! See above')
            sys.exit()
            
    def get_output_coords(self, iteration):        
        
        resultNames = ['Coord'] # 'Coord', 'U', 'RF'
        with open(self.workPath+'\\'+self.fileResultNames, "w") as f:
            for iResult in resultNames:
                f.write(iResult+'\n')
        
        self.GetOdbResults(iteration)
        # Get all nodes in boundary for global
        for jR in resultNames:
            self.nodesCoord_odb = np.load(self.workPath+'\\'+jR+self.fileOutputName+'.npy')
            os.remove(self.workPath+'\\'+jR+self.fileOutputName+'.npy')
            
        np.save(self.workPath+'Coord_odb.npy',self.nodesCoord_odb)
                
    def get_react_force_react_moment(self,initialRF, iteration):
        odbName = self.odbName[:-4]+'_job'+str(iteration)+'.odb'
        open(self.workPath+'\stderr.txt', 'w').close() #erase contents of error file
        with open(self.workPath+'\stdFlag.txt','ab') as out, open(self.workPath+'\stderr.txt','ab') as err:
            abaqusCall1 = 'cmd.exe /c '+self.abaqusVersion+' cae noGUI=GetGlobalRF.py' \
                            +' -- '+self.fileOutputName+' '+initialRF+' '+self.globLocDomain+' '+self.boundary+' ' \
                            +self.stepName+' '+odbName+' '+self.workPath
            process1 = subprocess.Popen(abaqusCall1, cwd=self.scriptPath,stdout=out,stderr=err)
            process1.wait()
            process1.terminate()
        if os.stat(self.workPath+'\stdFlag.txt').st_size != 0:
            with open(self.workPath+'\\'+'stdFlag.txt') as f:
                for cnt, line in enumerate(f):
                    print(line)
            with open(self.workPath+'\\'+'stderr.txt') as f:
                for cnt, line in enumerate(f):
                    print(line)
            self.printShell('Abaqus Error see above')
            sys.exit()
        
        # Get all nodes in boundary for global
        nodesRF_odb = -np.load(self.workPath+'\\'+self.fileOutputName+'RF'+'.npy')
        
        #self.printShell('Global reaction')
        
        os.remove(self.workPath+'\\'+self.fileOutputName+'RF'+'.npy')
        
        return nodesRF_odb
                
        
    def get_displacement_rotation(self,iteration):
        
        resultNames = ['UT'] # 'Coord', 'U', 'RF'
        with open(self.workPath+'\\'+self.fileResultNames, "w") as f:
            for iResult in resultNames:
                f.write(iResult+'\n')
        
        self.GetOdbResults(iteration)
        # Get all nodes in boundary for global
        nodesDispl_odb = np.load(self.workPath+'\\'+'UT'+self.fileOutputName+'.npy')
        os.remove(self.workPath+'\\'+'UT'+self.fileOutputName+'.npy')
        
        resultNames = ['UR'] # 'Coord', 'U', 'RF'
        with open(self.workPath+'\\'+self.fileResultNames, "w") as f:
            for iResult in resultNames:
                f.write(iResult+'\n')
        
        self.GetOdbResults(iteration)
        # Get all nodes in boundary for global
        nodesRot_odb = np.load(self.workPath+'\\'+'UR'+self.fileOutputName+'.npy')
        os.remove(self.workPath+'\\'+'UR'+self.fileOutputName+'.npy')
        
        #self.printShell('Global displ and rot')
        #self.printShell(nodesDispl_odb[:10])
        #self.printShell(nodesRot_odb[:10])
        
        return np.append(nodesDispl_odb,nodesRot_odb,axis=1)
        
    def get_odb_RF(self,iteration):
        
        resultNames = ['RF','RM'] # 'Coord', 'U', 'RF'
        with open(self.workPath+'\\'+self.fileResultNames, "w") as f:
            for iResult in resultNames:
                f.write(iResult+'\n')
        
        self.GetOdbResults(iteration)
        # Get all nodes in boundary for global
        nodesRF_odb = []
        for jR in resultNames:
            if jR == 'RF':
                nodesRF_odb = np.load(self.workPath+'\\'+jR+self.fileOutputName+'.npy')
                os.remove(self.workPath+'\\'+jR+self.fileOutputName+'.npy')
            elif jR == 'RM':
                nodesRF_odb = np.concatenate((nodesRF_odb,
                                              np.load(self.workPath+'\\'+jR+self.fileOutputName+'.npy')),
                                              axis = 1)
                os.remove(self.workPath+'\\'+jR+self.fileOutputName+'.npy')
        return nodesRF_odb
    
    def get_loc_SIF(self,iteration):
        
        t_odbName = self.odbName[:-4]+'_job'+str(iteration)+'.odb'
        open(self.workPath+'\\stderr.txt', 'w').close() #erase contents of error file
        with open(self.workPath+'\\stdFlag.txt','wb') as out, open(self.workPath+'\\stderr.txt','ab') as err:
            abaqusCall1 = 'cmd.exe /c '+self.abaqusVersion+' cae noGUI=GetOdbSIFResults.py' \
                            +' -- '+self.fileOutputName+' '+self.stepName+' '+t_odbName+' '+self.workPath
            process1 = subprocess.Popen(abaqusCall1, cwd=self.scriptPath,stdout=out,stderr=err)
            process1.wait()
            process1.terminate()
            del process1
        if os.stat(self.workPath+'\\stdFlag.txt').st_size != 0:
            with open(self.workPath+'\\'+'stdFlag.txt') as f:
                for cnt, line in enumerate(f):
                    print(line)
            with open(self.workPath+'\\'+'stderr.txt') as f:
                for cnt, line in enumerate(f):
                    print(line)
            self.printShell('Abaqus Error! See above')
            sys.exit()
        SIFResults = np.load(self.workPath+'\\SIF'+self.fileOutputName+'.npy')
        return SIFResults
    