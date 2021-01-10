#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 15:15:31 2020

@author: ch199899
"""

import os
import sys
import tempfile
import shutil
import numpy as np
from subprocess import call


## SCRIPT params
subjectName = 'CHOA_RENAL_006'         # name of the stud/center
mrn = 1000006                           # fake mrn for the study. Started it at 1000000 and keeps adding +1 ... in need for a better convention
scanDate = '20200803'                   # date of the scan. Not really important, but keeps consistency
folderTag = 'STARVIBE_dyn'     # folder tag. This indicates the characters that identify the folders with valid data. In theory not necessary when the DICOm folder ONLY contains the right folders, but just in case...
tempResol = 3.18                        # temporal resolution (to be used in the tags, not really important I think)
T = 156                                 # number of volumes (i.e. images)

## general path details
multicenterDicomFolder = '/fileserver/external/rdynaspro4/abd/multiCenterStudy/RAW/' + subjectName + '/DICOM/'
multicenterNRRDmFolder = '/fileserver/external/rdynaspro4/abd/multiCenterStudy/RAW/' + subjectName + '/NRRD/'
dataRawPath = '/fileserver/commondataraw/MRUoriginal/'
reconSCANPath = '/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/' + subjectName + '/'

## create new folder
dataRawPath += subjectName + '/' + str(mrn) + '/' +  scanDate + '/'
if not os.path.exists(dataRawPath):
    os.makedirs(dataRawPath)

if not os.path.exists(reconSCANPath):
    os.makedirs(reconSCANPath)

## convert DCM to NRRD
call(['crlDcmConvert','-fr',multicenterDicomFolder,multicenterNRRDmFolder ])

## copy files from NRRD to reconSCAN folder
files = os.listdir( multicenterNRRDmFolder )
ix = np.argsort(np.array([ int(ff.split('.')[0]) for ff in files ]))
joinCommand = ['/opt/el7/pkgs/crkit/release-current/bin/crlConvertN3DTo4D']
strTT = 0
for tt in range(T):
    
    if files[ix[tt]].find( folderTag ) > -1:
        subFiles = os.listdir(multicenterNRRDmFolder + files[ix[tt]])
        for ff in subFiles:
            if ff.endswith('.nrrd'):
                oldFile = multicenterNRRDmFolder + files[ix[tt]] + '/' + ff
    
        newFile = reconSCANPath + 'reconSCAN_T%d.nii' %(strTT)
        call(['crlConvertBetweenFileFormats', '-in',oldFile, '-out', newFile])
        
        strTT += 1
    
        joinCommand += ['-i', newFile]
        
joinCommand += ['-o',  reconSCANPath + 'reconSCAN_4D.nii']
call(joinCommand)
        
#  Include the info in the excel files manually and run the rest of the script 
#   runfile('/fileserver/external/rdynaspro4/abd/MRUcommon/singlePatientFullProcessing.py', wdir='/fileserver/external/rdynaspro4/abd/MRUcommon/')
        
     