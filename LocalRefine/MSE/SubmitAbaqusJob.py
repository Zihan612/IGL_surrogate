# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 09:32:04 2020

@author: RDCHLTBF
"""

from abaqus import mdb
import sys
import os

workPath = sys.argv[-2]
jobName = sys.argv[-1]

fPrint = open(workPath+'\\stderr.txt','a') 
sys.stdout=fPrint

path = workPath + '\\' + jobName
print(path)
os.chdir(workPath)

jobname = mdb.JobFromInputFile(jobName[:-4], path) 
jobname.submit() 
jobname.waitForCompletion()