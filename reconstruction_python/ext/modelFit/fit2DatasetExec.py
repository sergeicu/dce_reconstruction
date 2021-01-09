#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 09:13:29 2018

@author: jaume
"""
#%% IMPORTS
# general python imports
import os
import sys
import shutil
from subprocess import call
import argparse
import tempfile
from shutil import copyfile


# math importd
import numpy as np
import hdf5storage
import scipy.io as spio
from scipy.ndimage import morphology as mrp
import nibabel as nib
import nrrd


from rLTI_modelfit import runLTISysID4DCE


## MAIN FUNCTION
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='')
    parser.add_argument('-i', '--input', required=True,
                        help='input mat file to grasp data')
    parser.add_argument('-o', '--output', required=True,
                        help='output mat file')
    parser.add_argument('-mA', '--amask', required=True,
                        help='path to aorta mask')
    parser.add_argument('-mLK', '--lkmask', required=None,
                        help='path to left kidney mask')
    parser.add_argument('-mRK', '--rkmask', required=None,
                        help='path to right kidney mask')
    parser.add_argument('-strtFrame', '--endBaselineIX', default=0,
                        help='Starting frame of the injection')
    parser.add_argument('-tau', '--tau', default=100.0,
                        help='Regularization parameter Tau (default 100.0)')
    parser.add_argument('-numAtoms', '--numAtoms', default=100,
                        help='Number of atoms sampled in each iteration (default 100)')
    parser.add_argument('-numIter', '--numIter', default=100,
                        help='Number of inters to run (default 1000)')
    parser.add_argument('-numCluster', '--numCluster', default=100,
                        help='Number of clusters to use in each slice (default 100)')
    parser.add_argument('-p', '--cores', default=20,
                        help='Number of cores to use (default 20)')
    args = parser.parse_args()


    num_cores = int( args.cores )
    numAtoms = int( args.numAtoms )
    tau = int( args.tau )
    numClusters = int(args.numCluster)
    maxIter = int(args.numIter)
    endBaselineIX = int(args.endBaselineIX)

    dilation = np.ones([20,20,1])

    badVolumes = []

    #%% paths
    graspFile = args.input
    aortaMask = args.amask
    rightKMask = args.rkmask
    leftKMask = args.lkmask

    if not os.path.exists(os.path.dirname(args.output)):
        os.mkdir(os.path.dirname(args.output))

    #%% LOAD data
    print 'Loading: ' + graspFile
    graspData = np.abs(hdf5storage.loadmat(graspFile)['rec_image_GRASP'].transpose(0,1,3,2))
    

    #%% load masks
    masks = {'aorta':[]}
    #masks['aorta'] = nib.load(aortaMask).get_fdata()
    masks['aorta'] = nrrd.read(aortaMask)[0][:,::-1,:].transpose(0,2,1)
    masks['aifMask'] = nrrd.read(aortaMask)[0][:,::-1,:].transpose(0,2,1)

    fullMask = masks['aorta']
    if not rightKMask is None:
        #masks['right'] = nib.load(rightKMask).get_fdata()
        masks['right'] = nrrd.read(rightKMask)[0][:,::-1,:].transpose(0,2,1)
        masks['right'] = mrp.binary_dilation( masks['right'], dilation)
        fullMask += masks['right']

    if not leftKMask is None:
        #masks['left'] = nib.load(leftKMask).get_fdata()
        masks['left'] = nrrd.read(leftKMask)[0][:,::-1,:].transpose(0,2,1)
        masks['left'] = mrp.binary_dilation( masks['left'], dilation)
        fullMask += masks['left']
    

    ## create full masks
    masks['tissue'] = np.ones( masks['aorta'].shape )
    masks['tissue'] [np.where(fullMask)] = 0


    #%% run fitting
    ltiImg, aif = runLTISysID4DCE( graspData, masks, numAtoms=numAtoms, tau=tau, badVolumes=badVolumes, maxIter=maxIter, numClusters=numClusters, nthreads=num_cores, aif=None, endBaselineIX=endBaselineIX  )

    spio.savemat(args.output,{'ltifit':ltiImg, 'aif':aif})
