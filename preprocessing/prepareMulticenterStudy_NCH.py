
# copied from: /fileserver/external/rdynaspro4/abd/MRUcommon/jcollfont/ by sv407

import os
import sys
import tempfile
import shutil
from subprocess import call


## SCRIPT params
subjectName = 'NCH_PR_CMRU_001'         # name of the stud/center
mrn = 1000007                           # fake mrn for the study. Started it at 1000000 and keeps adding +1 ... in need for a better convention
scanDate = '20200220'                       # date of the scan. Not really important, but keeps consistency
folderTag = 'IM-0001'                        # folder tag. This indicates the characters that identify the folders with valid data. In theory not necessary when the DICOm folder ONLY contains the right folders, but just in case...
tempResol = 10.1                           # temporal resolution (to be used in the tags, not really important I think)
T = 70                                 # number of volumes
filetype = '.dcm'  # serge - added to account for non folder structure. If '.dcm' - we will loop over files, if 'folder' - we will loop over dicom folders 

## general path details
multicenterDicomFolder = '/fileserver/external/rdynaspro4/abd/multiCenterStudy/RAW/' + subjectName + '/DICOM/'
dataRawPath = '/fileserver/commondataraw/MRUoriginal/'
reconSCANPath = '/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/' + subjectName + '/'

## create new folder
dataRawPath += subjectName + '/' + str(mrn) + '/' +  scanDate + '/'
if not os.path.exists(dataRawPath):
    os.makedirs(dataRawPath)

## copy original DICOM files to commondataraw folder
folders = os.listdir(multicenterDicomFolder)
folders.sort()
newFolders = []
tt = 1                                                  # counter for sorting folders
for ff in folders:
    if ff.find(folderTag) > -1:
        newFolders.append( dataRawPath + str(tt) + '_COR_DYN_STARVIBE_dyn_TT=%0.2f/' %( tempResol*(tt-1) ) )
        shutil.copytree( multicenterDicomFolder + ff, newFolders[-1] )
        tt += 1                                         # update counter


## folder structure has been set. Now run the following python script without downloading data:
#        /fileserver/external/rdynaspro4/abd/MRUcommon/singlePatientFullProcessing.py
#        Remember to run it in a server with GPUs to compute the masks
#       Also, the script will likely fail when trying to get infor from the DICOM files (line 288):
#           'x=dicomHeaderInfoExtractor(patientName,patPath, seqPrefix1, seqPrefix2,numSeq);'
#        When it fails run the following:
# 
#   runfile('/fileserver/external/rdynaspro4/abd/MRUcommon/singlePatientFullProcessing.py', wdir='/fileserver/external/rdynaspro4/abd/MRUcommon/')

## resample all images to standard resolution and create 4D file and correct the orientation change that the resampling function applies (...)
joinCommand = ['/opt/el7/pkgs/crkit/release-current/bin/crlConvertN3DTo4D']
for tt in range(T):
    currentFile = reconSCANPath + 'reconSCAN_T%d.nii' %(tt)
    call([ 'crlResampler2', '--voxelsize', '1.25,3.0,1.25', '-i', currentFile, '-o', currentFile ])
    call([ 'crlOrientImage', currentFile, currentFile, 'coronal' ])
#    call([ 'crlCropImage', '-x', '32,0,0,224,224,36', '-i', currentFile, '-o', currentFile ])
    joinCommand += ['-i', currentFile]
joinCommand += ['-o',  reconSCANPath + 'reconSCAN_4D.nii']
call(joinCommand)

        
        
#  Include the info in the excel files manually and run the rest of the script 
#   runfile('/fileserver/external/rdynaspro4/abd/MRUcommon/singlePatientFullProcessing.py', wdir='/fileserver/external/rdynaspro4/abd/MRUcommon/')
        
        