#!/usr/bin/env python2# -*- coding: utf-8 -*-
"""

@author: jaume Coll-Font <jcollfont@gmail.com>
"""

# general python imports
import os
import shutil
from subprocess import call
import sys
import argparse
from joblib import Parallel, delayed
#import matplotlib.pyplot as plt

# math importd
import numpy as np
from scipy import ndimage 
import nrrd

import xlrd
import nibabel as nib


def parseXLSforSubjectData( patientName, xlsFile ):


    numPatients  = len(patientName)

    ## read subjects list file containing patient name/date and address and name of .dat files 
    workbook = xlrd.open_workbook(xlsFile)
    sheet = workbook.sheet_by_index(0)

    # if no patient specified, take them all!
    if numPatients == 0:
        patientName = sheet.row_values(0)[2:]

    subjectInfo = []
    for pt in patientName:

        subjectInfo.append( {'patientName':pt, 'MRN':-1,'files':[],'timePoints':[]} )


        colIX = sheet.row_values(0).index(pt)
        patientInfo = sheet.col_values(colIX)

        numOfDatFiles = 3 - int(patientInfo[3:6].count(''))
        
        subjectInfo[-1]['files'] = patientInfo[3:(3+numOfDatFiles)]           # retrieve .dat files with K-space data (possibly more than one) 
        subjectInfo[-1]['timePoints'] = patientInfo[7:(7+numOfDatFiles)]    # retrieve the number of time points / volumes to use

        subjectInfo[-1]['MRN'] = int(patientInfo[2])


    return  subjectInfo



#%%
class dciException(Exception):
    pass


## MAIN FUNCTION
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='FAST TLE algorithm for DCI estimation')
    parser.add_argument('-x', '--xls', required=True,
                        help='Excel file with all the subjects info')
    parser.add_argument('-s', '--subject', default='',
                        help='Subject Name')
    args = parser.parse_args()

    # prepare data
    patientsName = []
    if not args.subject == '':
        patientsName = [args.subject]

    # extract info from XLS
    subjectsInfo = parseXLSforSubjectData( patientsName, args.xls )

    # for all subjects
    for subj in subjectsInfo:

        # get paths to output filename and available segmentations
        coil_folder = 'data_DCE/DCE/Processed/Subj%d/%s/common_processed/coil_profile/' \
                                %( subj['MRN'], subj['patientName'] )
        
        virtual_coilweigths_filename = coil_folder + 'subj%d_virtual_coil_weigths.npz' \
                                            %( subj['MRN'] )
        virtual_coil_filename = coil_folder + 'subj%d_virtual_coil.nrrd' \
                                            %( subj['MRN'] )
        # try:   

        # load kidney masks
        left_kidney_mask_filename = '/fileserver/abd/GraspRecons/reconResultsSCAN/%s_seq1/leftKidneyMask.nii' \
                                        %(subj['patientName'])
        right_kidney_mask_filename = '/fileserver/abd/GraspRecons/reconResultsSCAN/%s_seq1/rightKidneyMask.nii' \
                                        %(subj['patientName'])
        aorta_mask_filename = '/fileserver/abd/GraspRecons/reconResultsSCAN/%s_seq1/aortaMask.nii' \
                                        %(subj['patientName'])
        
        
        left_kidney_mask = nib.load(left_kidney_mask_filename).get_fdata()
        right_kidney_mask = nib.load(right_kidney_mask_filename).get_fdata()
        aorta_mask = nib.load(aorta_mask_filename).get_fdata()

        # join all masks
        targetMask = left_kidney_mask + right_kidney_mask + aorta_mask
        targetMask[targetMask>0] = 1
        targetMask  = ndimage.zoom(targetMask, (2,2,19/float(targetMask.shape[2])) , order=0 )

        # load precomputed coil profile
        coilprofile = []
        for cc in os.listdir(coil_folder):
            if cc.find('nufft_coil_profile') > -1  & cc.endswith('.nrrd'):
                print cc
                coilprofile.append( nrrd.read( coil_folder + cc )[0] )
    
        coilprofile = np.array(coilprofile)

        # load 
        # RootSumSquaresCoil = np.sqrt(np.sum( coilprofile**2, axis=0 ))
    
        # solve LSQR problem
        S = coilprofile.reshape( coilprofile.shape[0], np.prod(coilprofile.shape[1:]) ).T
        b = targetMask#*RootSumSquaresCoil
        print S.shape
        print b.ravel().shape
        w = np.linalg.lstsq( S, b.ravel() )
        print w
        print len(w)
        print w[0].shape

        # save weights
        np.savez(virtual_coilweigths_filename, w)
        
        # compute virtual coil
        virtual_coil = np.zeros(coilprofile.shape[1:])
        for ii in range(coilprofile.shape[0]):
            virtual_coil = w[0][ii] * coilprofile[ii,:,:,:]
        
        nrrd.write(virtual_coil_filename,virtual_coil)

        # except dciException as e:
        #     errorPrefix = 'ERROR while running weigths for subject %d' %( subj['MRN'] )
        #     print( errorPrefix + str(e))
        #     exit(1)
