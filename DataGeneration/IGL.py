# -*- coding: utf-8 -*-
"""
Created on Mon Nov 30 12:54:31 2020

@author: RDCHLTBF
"""

import sys
import os
from ShellAbaqus import ShellAbaqus
import numpy as np
import subprocess
import time
import pandas as pd

class IGL:
    def __init__(self,info):
    #def __init__(self,workPath,scriptPath,info):
        
        #self.workPath = workPath
        #self.scriptPath = scriptPath
        self.index = info['index']
        self.abaqusVersion = info['IGLInfo']['abaqusVersion']
        self.locPath = info['locPath']
        self.scriptPath = info['scriptPath']
        IGLInfo = info['IGLInfo']
        locInfo = info['locInfo']
        globInfo = info['globInfo']
        self.globType = IGLInfo['globType']
        self.locType = IGLInfo['locType']
        self.algorithm = IGLInfo['algorithm']
        self.tol = IGLInfo['tol']
        self.maxSteps = IGLInfo['maxSteps']
        self.SIFConvergence = IGLInfo['SIFConvergence']
        
        locDir = os.listdir(self.locPath)
        for item in locDir: # Remove all lock files in directory
            if item.endswith(".lck"):
                os.remove(os.path.join(self.locPath, item))        
        
        
        t0 = time.time()
        globInfo['inpName'] = globInfo['CAEName'][:-4]+'.inp'
        self.Glob = ShellAbaqus(self.abaqusVersion,self.locPath,self.scriptPath,globInfo)

        t1 = time.time()
        self.print_IGL('Global initiation time is '+str(t1-t0))
        locInfo['inpName'] = locInfo['CAEName'][:-4]+'.inp'
        self.Loc = ShellAbaqus(self.abaqusVersion,self.locPath,self.scriptPath,locInfo)
        
        t2 = time.time()
        print('Local initiation time is '+str(t2-t1))
        
        #
        self.abaqusVersion=IGLInfo['abaqusVersion']
        self.globInpName = globInfo['CAEName'][:-4]+'.inp'
        self.locInpName = locInfo['CAEName'][:-4]+'.inp'
        
        self.LocTime = 0.
        self.GlobTime = 0.
        self.IGLTime = 0.
        
        self.iteration = 0
        self.error = [10000.]
        if self.locType == 'AbaqusSubStructure' and self.SIFConvergence == True:
            self.print_IGL('Error: Cannot finds SIFs with substructured local problem.Please change SIFConvergence to False.')
            sys.exit()
        self.LocalSIF = {}
        self.GlobalExtraLoad = {} # Dictionaries with iterations as keys
        self.GlobalRF = {}
        self.GlobalTraction = {} 
        self.LocalAppliedDispl = {}
        self.LocalReactionTrain = {}
        self.LocalReaction = {}
        self.AuxReaction = {}
        self.Residual = {}
        self.mapGlob2Loc = []
        self.mapLoc2Glob = []
        self.mapGlob2Aux = []
        self.mapAux2Glob = []
        self.mapGlob2Glob = []
        self.mapGlob2Glob_2 = []
        open(self.locPath+'\stdout.txt', 'w').close() #erase contents of output file
    def print_IGL(self,text):
        with open(self.locPath+'stdout.txt', 'a') as f:
            print(text, file=f)
    
    def mapping_global_applied_displ_2_local_nodes(self):
        # ind = np.lexsort((first_names, surnames))
        
        globInd = list(range(self.Glob.nodesCoord_odb.shape[0]))
        self.mapGlob2Loc = []
        for i in self.Loc.nodesCoord_inp: #self.Glob.nodesCoord_odb:
            for j,jInd in enumerate(globInd):
                distance = np.linalg.norm(self.Glob.nodesCoord_odb[jInd] - i)
                if distance < 10e-4:
                    self.mapGlob2Loc += [jInd]
                    del globInd[j]
                    break
                if j == (len(globInd)-1):
                    self.print_IGL(str(self.Glob.nodesCoord_odb))
                    self.print_IGL('Something wrong with global model. Non-matching mesh: no node found for '+str(i))
                    sys.exit()
        if self.locType == 'GaussianProcess':
            mapLoc2Glob = []
            locInd = list(range(self.Loc.nodesCoord_inp.shape[0]))
            for i in self.Glob.nodesCoord_odb: #self.Glob.nodesCoord_odb:
                for j,jInd in enumerate(locInd):
                    distance = np.linalg.norm(self.Loc.nodesCoord_inp[jInd] - i)
                    if distance < 10e-4:
                        mapLoc2Glob += [jInd]
                        del locInd[j]
                        break
                    if j == (len(locInd)-1):
                        self.print_IGL(str(self.Loc.nodesCoord_inp))
                        self.print_IGL('Something wrong with local model. Non-matching mesh: no node found for '+str(i))
                        sys.exit()
            self.Loc.mapLoc2Glob = mapLoc2Glob
    
    def mapping_local_reaction_2_global_inp_nodes(self):
        
        self.mapLoc2Glob = []
        locInd = list(range(self.Loc.nodesCoord_odb.shape[0]))
        for i in self.Glob.nodesCoord_inp:
            for j,jInd in enumerate(locInd):
                distance = np.linalg.norm(self.Loc.nodesCoord_odb[jInd] - i)
                if distance < 10e-4:
                    self.mapLoc2Glob += [jInd]
                    del locInd[j]
                    break
                if j == (len(locInd)-1):
                    self.print_IGL(str(self.Loc.nodesCoord_odb))
                    self.print_IGL('Non-matching mesh: no node found for '+str(i))
                    sys.exit()
                    
        if self.locType == 'GaussianProcess':
            mapGlob2Loc = []
            globInd = list(range(self.Glob.nodesCoord_inp.shape[0]))
            for i in self.Loc.nodesCoord_odb:
                for j,jInd in enumerate(globInd):
                    distance = np.linalg.norm(self.Glob.nodesCoord_inp[jInd] - i)
                    if distance < 10e-4:
                        mapGlob2Loc += [jInd]
                        del globInd[j]
                        break
                    if j == (len(globInd)-1):
                        self.print_IGL(str(self.Glob.nodesCoord_inp))
                        self.print_IGL('Non-matching mesh: no node found for '+str(i))
                        sys.exit()
            self.Loc.mapGlob2Loc = mapGlob2Loc
                    
    def mapping_global_applied_displ_2_aux_nodes(self):
        # ind = np.lexsort((first_names, surnames))
        
        globInd = list(range(self.Glob.nodesCoord_odb.shape[0]))
        self.mapGlob2Aux = []
        for i in self.Aux.nodesCoord_inp: #self.Glob.nodesCoord_odb:
            for j,jInd in enumerate(globInd):
                distance = np.linalg.norm(self.Glob.nodesCoord_odb[jInd] - i)
                if distance < 10e-4:
                    self.mapGlob2Aux += [jInd]
                    del globInd[j]
                    break
                if j == (len(globInd)-1):
                    self.print_IGL(str(self.Glob.nodesCoord_odb))
                    self.print_IGL('Non-matching mesh: no node found for '+str(i))
                    sys.exit()
    
    def mapping_aux_reaction_2_global_inp_nodes(self):
        
        self.mapAux2Glob = []
        locInd = list(range(self.Aux.nodesCoord_odb.shape[0]))
        for i in self.Glob.nodesCoord_inp:
            for j,jInd in enumerate(locInd):
                distance = np.linalg.norm(self.Aux.nodesCoord_odb[jInd] - i)
                if distance < 10e-4:
                    self.mapAux2Glob += [jInd]
                    del locInd[j]
                    break
                if j == (len(locInd)-1):
                    self.print_IGL(str(self.Aux.nodesCoord_odb))
                    self.print_IGL('Non-matching mesh: no node found for '+str(i))
                    sys.exit()
    
    def mapping_global_traction_2_global_inp_nodes(self):
        
        self.mapGlob2Glob = []
        #for key, value in self.nodesCoord_inp.items():
        globInd = list(range(self.Glob.nodesCoord_odb.shape[0]))
        for i in self.Glob.nodesCoord_inp:
            for j,jInd in enumerate(globInd):
                distance = np.linalg.norm(self.Glob.nodesCoord_odb[jInd] - i)
                if distance < 10e-4:
                    self.mapGlob2Glob += [jInd]
                    del globInd[j]
                    break
                if j == (len(globInd)-1):
                    self.print_IGL(str(self.Glob.nodesCoord_odb))
                    self.print_IGL('Non-matching mesh: no node found for '+str(i))
                    sys.exit() 

    def mapping_global_inp_nodes_2_global_traction(self):
        
        self.mapGlob2Glob_2 = []
        #for key, value in self.nodesCoord_inp.items():
        globInd = list(range(self.Glob.nodesCoord_inp.shape[0]))
        for i in self.Glob.nodesCoord_odb:
            for j,jInd in enumerate(globInd):
                distance = np.linalg.norm(self.Glob.nodesCoord_inp[jInd] - i)
                if distance < 10e-4:
                    self.mapGlob2Glob_2 += [jInd]
                    del globInd[j]
                    break
                if j == (len(globInd)-1):
                    self.print_IGL(str(self.Glob.nodesCoord_inp))
                    self.print_IGL('Non-matching mesh: no node found for '+str(i))
                    sys.exit() 
        
    def fixed_point(self):
        # This assumes an .inp file is in work directory
        self.GlobalExtraLoad[self.iteration] = np.zeros((len(self.Glob.nodesInd_inp),6))
        t0 = time.time()
        while self.error[-1] > self.tol:
            # Solve global model with extra load, extract displacement and traction
            self.Glob.apply_extra_load(self.GlobalExtraLoad[self.iteration],self.iteration)
            a = time.time()
            self.Glob.run_job(self.iteration)
            b = time.time()
            self.GlobTime = self.GlobTime + b - a
            if self.iteration == 0:
                self.Glob.get_output_coords(self.iteration)
                self.mapping_global_traction_2_global_inp_nodes()
                self.mapping_global_applied_displ_2_local_nodes()
                self.GlobalRF[self.iteration] = self.Glob.get_react_force_react_moment('True',self.iteration)[self.mapGlob2Glob]
            else:
                self.GlobalRF[self.iteration] = self.Glob.get_react_force_react_moment('False',self.iteration)[self.mapGlob2Glob]
            self.LocalAppliedDispl[self.iteration] = self.Glob.get_displacement_rotation(self.iteration)[self.mapGlob2Loc]
            self.GlobalTraction[self.iteration] = self.GlobalRF[self.iteration] - self.GlobalExtraLoad[self.iteration]
            
            # Solve Local model with imposed displacement, extract reaction
            self.Loc.apply_displacement(self.LocalAppliedDispl[self.iteration],self.iteration)
            a = time.time()
            self.Loc.run_job(self.iteration)
            b = time.time()
            self.LocTime = self.LocTime + b - a
            if self.iteration == 0:
                self.Loc.get_output_coords(self.iteration)
                self.mapping_local_reaction_2_global_inp_nodes()
            self.LocalReaction[self.iteration] = self.Loc.get_odb_RF(self.iteration)[self.mapLoc2Glob]
            
                
            # Compute residual
            self.Residual[self.iteration] = (self.GlobalTraction[self.iteration] - self.LocalReaction[self.iteration])
            
            # Update global extra load
            #self.GlobalExtraLoad[self.iteration+1] = self.GlobalExtraLoad[self.iteration] + self.Residual[self.iteration]
            self.GlobalExtraLoad[self.iteration+1] = self.LocalReaction[self.iteration]
            
            if self.SIFConvergence == True and self.locType == 'Abaqus':
                self.LocalSIF[self.iteration] = self.Loc.get_loc_SIF(self.iteration)
                self.error += [np.max(np.abs(self.Residual[self.iteration]))]
                '''
                if self.iteration == 0:
                    self.error += [self.error[0]]
                else:
                    self.error += [np.max(np.abs((self.LocalSIF[self.iteration]-self.LocalSIF[self.iteration-1])/self.LocalSIF[self.iteration]))]
                '''
            else:
                self.error += [np.max(np.abs(self.Residual[self.iteration]))]
            self.print_IGL('Iteration = '+str(self.iteration))
            self.print_IGL(str(self.error))
            self.iteration += 1
            
            if self.iteration > self.maxSteps:
                self.print_IGL('IGL did not converge.')
                sys.exit()
        
        t1 = time.time()
        self.IGLTime = self.IGLTime + t1 - t0
            
    def fixed_point_aitken(self):
        self.print_IGL('Started Aitken IGL')
        # This assumes an .inp file is in work directory
        Aitken = 1.
        self.GlobalExtraLoad[self.iteration] = np.zeros((len(self.Glob.nodesInd_inp),6))
        t0 = time.time()
        while self.error[-1] > self.tol:
            # Solve global model with extra load, extract displacement and traction
            a = time.time()
            self.print_IGL('Global apply extra load')
            self.Glob.apply_extra_load(self.GlobalExtraLoad[self.iteration],self.iteration)
            self.Glob.run_job(self.iteration)
            if self.iteration == 0:
                t01 = time.time()
                self.Glob.get_output_coords(self.iteration)
                self.mapping_global_traction_2_global_inp_nodes()
                self.mapping_global_applied_displ_2_local_nodes()
                t11 = time.time()
                self.print_IGL('Global coordinates and mapping takes '+str(t11-t01))
                self.GlobalRF[self.iteration] = self.Glob.get_react_force_react_moment('True',self.iteration)[self.mapGlob2Glob]
            else: 
                self.GlobalRF[self.iteration] = self.Glob.get_react_force_react_moment('False',self.iteration)[self.mapGlob2Glob]
            self.LocalAppliedDispl[self.iteration] = self.Glob.get_displacement_rotation(self.iteration)[self.mapGlob2Loc] # Saves LocalAppliedDispl[self.iteration]
            # df_u = pd.DataFrame(self.LocalAppliedDispl[self.iteration])
            # df_u.to_csv(self.scriptPath + '\\Files\OutputData\\LocalAppliedDispl\\S{}_Displ{}.csv'.format(self.index,self.iteration))
            b = time.time()
            self.GlobTime = self.GlobTime + b - a
            
            self.GlobalTraction[self.iteration] = self.GlobalRF[self.iteration] - self.GlobalExtraLoad[self.iteration]
            
            # Solve Local model with imposed displacement, extract reaction
            a = time.time()
            self.print_IGL('Applying local displacement')
            self.Loc.apply_displacement(self.LocalAppliedDispl[self.iteration],self.iteration)
            self.Loc.run_job(self.iteration)
            if self.iteration == 0:
                t21 = time.time()
                self.Loc.get_output_coords(self.iteration)
                np.save(self.locPath+'locCoordOdb',self.Loc.nodesCoord_odb)
                np.save(self.locPath+'locCoordInp',self.Loc.nodesCoord_inp)
                self.mapping_local_reaction_2_global_inp_nodes()
                t31 = time.time()
                self.print_IGL('Local coordinates and mapping takes '+str(t31-t21))
            self.LocalReactionTrain[self.iteration] = self.Loc.get_odb_RF(self.iteration) # Extra variable to allow training with surrogate models with mapping problems
            # df_p = pd.DataFrame(self.LocalReactionTrain[self.iteration])
            # df_p.to_csv(self.scriptPath + '\\Files\OutputData\LocalReaction\S{}_Reaction{}.csv'.format(self.index,self.iteration))
            self.LocalReaction[self.iteration] = self.LocalReactionTrain[self.iteration][self.mapLoc2Glob]
            b = time.time()
            self.LocTime = self.LocTime + b - a
                
            # Compute residual
            self.Residual[self.iteration] = self.GlobalTraction[self.iteration] - self.LocalReaction[self.iteration]
            
            # Update global extra load
            self.GlobalExtraLoad[self.iteration+1] = self.GlobalExtraLoad[self.iteration] + self.Residual[self.iteration] # Outputs Global Extra Load
            
            # Update global extra load
            if self.iteration > 0:
                tmp = (self.Residual[self.iteration]-self.Residual[self.iteration-1])
                numerator = np.sum(np.dot(self.Residual[self.iteration-1].T,tmp))
                denominator = np.sum(np.dot(tmp.T,tmp))
                Aitken = -Aitken * numerator/denominator
                
            self.GlobalExtraLoad[self.iteration+1] = self.GlobalExtraLoad[self.iteration] + Aitken * self.Residual[self.iteration]
            
            if self.SIFConvergence == True and self.locType == 'Abaqus':
                self.LocalSIF[self.iteration] = self.Loc.get_loc_SIF(self.iteration)
                self.error += [np.max(np.abs(self.Residual[self.iteration]))]
                
                #if self.iteration == 0:
                #    self.error += [1]
                #else:
                #    self.error += [np.max(np.abs((self.LocalSIF[self.iteration]-self.LocalSIF[self.iteration-1])/self.LocalSIF[self.iteration]))]
            else:
                self.error += [np.max(np.abs(self.Residual[self.iteration]))]
            self.print_IGL('Iteration = '+str(self.iteration)+', Error='+str(self.error[-1]))
            self.iteration += 1
            if self.iteration > self.maxSteps:
                self.print_IGL('IGL did not converge. Problem with convergence!!!')
                break
        t1 = time.time()
        self.IGLTime = self.IGLTime + t1 - t0
            
    def apply_displacement(self):
        jobName = self.locInpName[:-4]+'_postprocess.inp'
        # self.Loc.get_nodes_input_NSC()
        # self.Loc.get_nodes_coord_input_NSC()
        # self.Loc.map_inp_2_NSCinp()
        with open(self.locPath+'\\'+self.locInpName) as fin, open(self.locPath+'\\'+jobName,'w') as fout:     
            BCCounter = 5
            for line in fin:
                fout.write(line)
                if line == '*End Instance\n':  # Create all of the nodesets for the boundary conditions
                    for NodeNum in self.Loc.nodesInd_inp:
                        fout.write('*Nset, nset=GL-'+str(NodeNum)+', instance='+self.Loc.instanceNames+'\n')
                        fout.write(str(NodeNum)+',\n')
                
                if line == '*Step, name='+ self.Loc.stepName +', nlgeom=NO\n': # Create boundary conditions after step definition
                    BCCounter = 0
                BCCounter = BCCounter + 1
                if BCCounter == 4:
                    fout.write('** BOUNDARY CONDITIONS\n')
                    for n,NodeNum in enumerate(self.Loc.nodesInd_inp):
                        fout.write('*Boundary\n')
                        fout.write('GL-'+str(NodeNum)+', 1, 1, '+ str(self.LocalAppliedDispl[self.iteration][n,0]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 2, 2, '+ str(self.LocalAppliedDispl[self.iteration][n,1]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 3, 3, '+ str(self.LocalAppliedDispl[self.iteration][n,2]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 4, 4, '+ str(self.LocalAppliedDispl[self.iteration][n,3]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 5, 5, '+ str(self.LocalAppliedDispl[self.iteration][n,4]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 6, 6, '+ str(self.LocalAppliedDispl[self.iteration][n,5]) +'\n')
        fout.close()
        
    def apply_extra_load(self):
        
    ##############################################
    # Apply global extra load values to .inp file
    ############################################## 
        jobName = self.globInpName[:-4]+'_postprocess.inp'
        # self.Glob.get_nodes_input_NSC()
        # self.Glob.get_nodes_coord_input_NSC()
        # self.Glob.map_inp_2_NSCinp()
        with open(self.locPath+'\\'+self.globInpName) as fin, open(self.locPath+'\\'+jobName,'w') as fout:     
            for line in fin:
                fout.write(line)
                if line == '*End Instance\n':  # Create all of the nodesets for the boundary conditions
                    for n, NodeNum in enumerate(self.Glob.nodesInd_inp):
                        NodeNum = self.Glob.nodesInd_inp2[self.Glob.mapNSCInp2Inp[n]]
                        fout.write('*Nset, nset=GL-'+str(NodeNum)+', instance='+self.Glob.instanceNames+'\n')
                        fout.write(str(NodeNum)+',\n')
                
                if line == '** LOADS\n':
                    for n, NodeNum in enumerate(self.Glob.nodesInd_inp):
                        NodeNum = self.Glob.nodesInd_inp2[self.Glob.mapNSCInp2Inp[n]]
                        fout.write('*Cload\n')
                        fout.write('GL-'+str(NodeNum)+', 1, '+ str(self.GlobalExtraLoad[self.iteration][n,0]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 2, '+ str(self.GlobalExtraLoad[self.iteration][n,1]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 3, '+ str(self.GlobalExtraLoad[self.iteration][n,2]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 4, '+ str(self.GlobalExtraLoad[self.iteration][n,3]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 5, '+ str(self.GlobalExtraLoad[self.iteration][n,4]) +'\n')
                        fout.write('GL-'+str(NodeNum)+', 6, '+ str(self.GlobalExtraLoad[self.iteration][n,5]) +'\n')
        fout.close()
        
    def run_job(self,jobName):
                    
        open(self.locPath+'\stderr.txt', 'w').close() #erase contents of error file
        with open(self.locPath+'\stdFlag.txt','wb') as out, open(self.locPath+'\stderr.txt','ab') as err:
            abaqusCall1 = 'cmd.exe /c '+self.abaqusVersion+' cae noGUI='+self.scriptPath+'\\SubmitAbaqusJob.py' \
                            +' -- '+self.locPath+' '+jobName
            process1 = subprocess.Popen(abaqusCall1, cwd=self.locPath,stdout=out,stderr=err)
            process1.wait()
            process1.terminate()
            del process1
            
        if os.stat(self.locPath+'\stdFlag.txt').st_size != 0:
            with open(self.locPath+'\\'+'stdFlag.txt') as f:
                for cnt, line in enumerate(f):
                    print(line)
            with open(self.locPath+'\\'+'stderr.txt') as f:
                for cnt, line in enumerate(f):
                    print(line)
            self.print_IGL('Abaqus Error! See above')
            sys.exit()
            
    def postprocessing(self):
        t0 = time.time()
        self.iteration = self.iteration-1
        self.apply_displacement()
        self.run_job(self.locInpName[:-4]+'_postprocess.inp')
        t1 = time.time()
        self.print_IGL('Local Postprocess time is '+str(t1-t0))
       
    def run(self):
        
        if self.algorithm == 'FixedPoint':
            self.fixed_point()
        elif self.algorithm == 'FixedPointAitken':
            self.fixed_point_aitken()
            
        # Postprocess
        self.postprocessing()
        
            
