#!/usr/bin/env python2# -*- coding: utf-8 -*-

"""
@author: jaume Coll-Font <jcollfont@gmail.com>
"""

# print function workaround 
from __future__ import print_function


# general python imports
import os
import sys
import shutil
from subprocess import call
import argparse
from shutil import copyfile
import tempfile

# math imports
import numpy as np
import nibabel as nib
import nrrd
import scipy.io as spio
import hdf5storage
#import matplotlib.pyplot as plt
from skimage.morphology import closing, opening

# parallel processing imports
from joblib import Parallel, delayed
import multiprocessing




# local imports
sys.path.append('ext/modelFit/') # copied from - sys.path.append('/home/ch199899/Documents/Research/DCE_all/modelFit/')
from resizeTools import increaseDimVolume,increaseDimMask

sys.path.append('ext/image_reconstruction/') # 
# copied from - sys.path.append('/home/ch199899/Documents/Research/DCE_all/image_reconstruction/') 
# NB - modified path of this append by removing '_XDgrasp' from the suffix of last dir 
from conjugateGradientRecon import conjugateGradientRecon

# Note that the registration2 imports above are NOT technically required since no bulk motion detection flag is used on input (aka we use `-no_bmd` flag in the input args, which is equivalent to running the matlab GRASP reconstruction)
#sys.path.append('ext/registration2/') # copied from - sys.path.append('/home/ch199899/Documents/Research/DCE_all/registration2/')
#from RLTIegistration import runRLTIregistration, find_reference_img, create_movement_masks
#from sequential_registration import computeInverseRegTransform, applyRegistration, register_images2reference

#sys.path.append('ext/pynufft') 



#%%
#### --------------------------  RUN SCRIPT  -------------------------- ####

#%%
#
#
#
#
def detectCorruptedSpokes(FIDsignal, spk_per_vol, startVol):

    #%% load FID data
    FIDs_all = spio.loadmat(FIDsignal)
    corrFID = np.max(1-FIDs_all['FID_corr'],axis=0)

    NSpk =  corrFID.size

    ## THRESHOLD THE FID AND CREATE MASKS
    if not 'corrFID_thr' in locals():
        # corrFID_thr =  min( 0.005,max( 0.1, np.std(corrFID[ int(spk_per_vol*startVol): ])*2.  )  )
        corrFID_thr = max( 0.005, min(0.1, np.median(corrFID[ int(spk_per_vol*startVol): ]) + np.std(corrFID[ int(spk_per_vol*startVol): ])))
        # corrFID_thr =  np.std(corrFID[ spk_per_vol*startVol: ])*2.
    movementIX = corrFID > corrFID_thr
    OutlierMask = np.zeros([NSpk])
    OutlierMask[movementIX.ravel()] = 1

    # filter
    OutlierMask_filt = opening( OutlierMask.reshape(OutlierMask.size,1), np.ones([1,15]) ).ravel()

    # eliminate anything before startVol images
    OutlierMask_filt[ :int(spk_per_vol*startVol) ] = 0

    # outlierMask_lowres = np.zeros([T])
    # for rr in range(T):
    #     aa = range(spk_per_vol*rr,spk_per_vol*(rr+1))
    #     if np.sum( OutlierMask_filt[aa] ) > 0:
    #         outlierMask_lowres[rr] = 1
        
    # rem_spokes = NSpk - T * spk_per_vol
    # vol_per_spk = 1 / float(spk_per_vol)    
    # fid_time = np.arange(  1 - vol_per_spk * spk_per_vol/2. , \
    #                                 T  + vol_per_spk * (spk_per_vol/2. +  rem_spokes ), \
    #                                 vol_per_spk )
    # mean_img_time = np.arange(T) +1
    # endBaselineFIDIX = endBaselineIX*spk_per_vol

    # outlierMask_lowres = np.zeros([T])
    # for rr in range(T):
    #     aa = range(spk_per_vol*rr,spk_per_vol*(rr+1))
    #     if np.sum( OutlierMask_filt[aa] ) > 0:
    #         outlierMask_lowres[rr] = 1

    return np.where( OutlierMask_filt.ravel() )[0], corrFID_thr, corrFID


#%%
#
#
#
def invertMasks( maskFiles, outFolder, sliceOverSampling, newImgSize=None, referenceNii=None, rescaleBool=False, flipDim=-1 ):

    inverseMaskFile = []
    for mm in maskFiles:

        # correct for possible mask names (basically for right kidney... sort of hard coded, I know)
        if os.path.exists(mm):
            # if os.path.exists(mm + '.gz'):
            #     mm += '.gz'
            # elif os.path.exists( mm[:-4] + 'UpperLower.nii'):
            #     mm = mm[:-4] + 'UpperLower.nii'
            # else:
            #     filesList = glob.glob( mm[:-4] + '*' )
            #     if len(filesList) > 0:
            #        mm = filesList[0]

            maskNii = nib.load(mm)

            inverseMaskFileName = outFolder + os.path.basename( mm.split('.')[0] + '_inverted' )

            maskImg = maskNii.get_fdata()
            inverseMaskFile.append( increaseDimMask( maskImg, inverseMaskFileName, referenceNii,  flipDim=flipDim, sliceOverSampling=sliceOverSampling, newImgSize=newImgSize, rescaleBool=rescaleBool )[1] )

    return inverseMaskFile

#%%
#
#
#
def computeInverseLTIregistration( outFolder, registrationTFM, ltiFitImages, inverseLTIVolumeFile , sliceOverSampling, newImgSize, num_cores=20):

    T = len(ltiFitImages)

    # define names
    inverseTransforms = []
    inverseLTIFitImages = []
    for tt in range(T):
        inverseTransforms.append(  outFolder + os.path.basename(registrationTFM[tt])[:-4] + '_inv.txt'  )
        inverseLTIFitImages.append( outFolder + ''.join(os.path.basename(ltiFitImages[tt]).split('.')[:-1]) + '_invRegistered_' + str(tt) +  '.nii.gz' )

    # compute inverse transform
    Parallel(n_jobs=num_cores)(delayed(computeInverseRegTransform) \
                (ltiFitImages[tt], registrationTFM[tt], inverseTransforms[tt], tt) for tt in range( T ) )

    # apply inverse transform
    Parallel(n_jobs=num_cores)(delayed(applyRegistration) \
                (ltiFitImages[tt], ltiFitImages[tt], inverseTransforms[tt], inverseLTIFitImages[tt], tt,  mode='elastix_dense') for tt in range( T ) )

    
    # flip Z dimension
    flippedLTIFile = invertMasks( inverseLTIFitImages, outFolder,   \
                                    sliceOverSampling=sliceOverSampling, newImgSize=newImgSize, referenceNii=ltiFitImages[0], rescaleBool=False, flipDim=1)
    
#    # save registered LTI fit as .mat
#    for tt in range(T):
#        nibData = nib.load( inverseLTIFitImages[tt] )
#        if tt == 0:
#            ltiVolume = np.zeros( nibData.get_fdata().shape + (T,))
#        ltiVolume[:,:,:,tt] =  nibData.get_fdata()
#        
#    hdf5storage.savemat( inverseLTIVolumeFile, ltiVolume.transpose(0,1,3,2) )   # saving as Nx,Ny,T,Nz
#
#
    return flippedLTIFile, inverseTransforms

##
#
#
#
#
def computeEndBaseline(rec_image_NUFFT, masks=None):

    T = len(rec_image_NUFFT)
    imgDims = rec_image_NUFFT[0].shape
    
    ## apply masks
    if not masks is None:
        fullMask = np.zeros( imgDims )
        for mm in range(len(masks)):
            maskNii = nib.load(masks[mm]).get_fdata()
            fullMask += maskNii
        fullMask[fullMask>0] = 1

        maskedImages = np.array(rec_image_NUFFT) * np.tile(fullMask,[T,1,1,1])
        
    else:
        maskedImages = np.array(rec_image_NUFFT)

    ## get average signal throughout the tissue
    meanAllTissue = np.mean( np.abs(maskedImages).reshape( T, np.prod(imgDims) ), axis=1).ravel()

    ## compute the first differences of the tissue
    allTissueDiff = np.diff(meanAllTissue)

    ## get the baseline section.
    #   Here we assume that the differential is peaked and this peak drives
    #   the mean. Hence. we can get the points until the signal is
    #   surpassed.
    threshold = np.median(allTissueDiff) + np.std(allTissueDiff)
    
    return np.where((allTissueDiff > threshold).ravel())[0][0] 

def find_reference_scan(args,patientName):
    
    # find a reference from the scanner (i.e. image reconstructed ON the siemens scanner)


    # if not specified - attempt to find the reference inside the CRL filesystem
    if args.referenceSCAN is None:
        referenceFileSCAN = '/fileserver/abd/GraspRecons/reconResultsSCAN/'+ patientName +'/reconSCAN_T0.nii'
        if not os.path.exists(referenceFileSCAN):
            referenceFileSCAN = '/fileserver/abd/GraspRecons/reconResultsSCAN/'+ patientName +'_seq1/reconSCAN_T0.nii'
            if not os.path.exists(referenceFileSCAN):
                referenceFileSCAN = '/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/'+ patientName +'/reconSCAN_T0.nii'
                if not os.path.exists(referenceFileSCAN):
                    referenceFileSCAN = '/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/'+ patientName +'_seq1/reconSCAN_T0.nii'
                else:
                    referenceFileSCAN = 'FILE NOT FOUND'
                    print('Reference scan not found!')
    else:
        # if reference is specified 
        
        referenceFileSCAN = args.referenceSCAN # if full path to reference file is specified          
        if not os.path.isfile(referenceFileSCAN):
            # if only a folder is specified and the patient has a single sequence 
            referenceFileSCAN = args.referenceSCAN + patientName + 'reconSCAN_T0.nii' # if single sequence 
            if not os.path.isfile(referenceFileSCAN): 
                # if only a folder is specified and the patient has a dual sequence 
                referenceFileSCAN = args.referenceSCAN + patientName + '_seq1/reconSCAN_T0.nii' # if multi sequence 
    return referenceFileSCAN
    
#%%
#### --------------------------  RUN SCRIPT  -------------------------- ####


## MAIN FUNCTION
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='')
    parser.add_argument('-d', '--baseFolder', required=True,
                        help='Base folder.')
    parser.add_argument('-s', '--mrn', required=True,
                        help='subject s MRN')
    parser.add_argument('-o', '--output_folder', required=True,
                        help='output folder')
    parser.add_argument('-r', '--referenceSCAN', default=None,
                        help='reference nii file from SCANNER recon')
    parser.add_argument('-m', '--masksPath', default='',
                        help='Path to masks folder. Should contain an aortaMask.nii as well as right/leftKidneyMask.nii')
    parser.add_argument('-l', '--lamReg', default=0.0001,
                        help='Regularization Lambda')
    parser.add_argument('-llti', '--lamRegLTI', default=None,
                        help='Regularization Lambda')
    parser.add_argument('-stkBinWidth', '--stkBinWidth', default='0.1',
                        help='Width of each bin in the stokes phase. \
                            If multiple phases are defined, separate by a comma. \
                            First corresponds to time, second to respiration phase.')
    parser.add_argument('-pn', '--patientName', default=None,
                        help='Patient name (used to create folder structure)')
    parser.add_argument('-strtFrame', '--endBaselineIX', default=None,
                        help='Starting frame of the injection')
    parser.add_argument('-sl', '--slices', default=None,
                        help='Selected slices. Show comma separated')
    parser.add_argument('-f', '--force', action='store_true',
                        help='force recomputation of all results')
    parser.add_argument('-no_bmd', '--bulkMotionDetector', action='store_false',
                        help='bulk motion detection activatied? Default True, otherwise all spokes will be used for reconstruction.')
    parser.add_argument('-R', '--regMethod', default=None,
                        help='Which in-the-loop registration method to use: None, \'standard\', \'sequential\'. ')
    parser.add_argument('-nufft', '--nufftBool', action='store_true',
                        help='Compute an unregularized solution to start.')
    parser.add_argument('-pyramInit', '--pyramInit', action='store_true',
                        help='Initialize algorithm with a hyerarchical optimization over time.')
    parser.add_argument('-sliceOVer', '--slcieOver', default=0.2,
                        help='Slice oversampling')
    parser.add_argument('-p', '--cores', default=20,
                        help='Number of cores to use')
    args = parser.parse_args()


    ## Load Subject MRN
    mrn = args.mrn  

    # create folder structure
    baseFolder = args.baseFolder
    baseFolder += '/' + 'Subj' + mrn + '/'
    if not os.path.exists(baseFolder):
        os.makedirs(baseFolder)
    if args.patientName is None:
        patientName = os.listdir(baseFolder)[0] 
    else:
        patientName = args.patientName
    baseFolder += patientName + '/'
    if not os.path.exists(baseFolder):
        os.makedirs(baseFolder)
    baseFolder += 'common_processed/'
    if not os.path.exists(baseFolder):
        os.makedirs(baseFolder)

    print('Running RLTI processing for subject: %s, MRN: %s' %( patientName, mrn ))
    referenceFileSCAN = find_reference_scan(args, patientName)
    
    # recon params
    num_cores = int(args.cores)
    lambdaBase = float( args.lamReg )
    badVolumes = []
    num_stokes_per_image = 34
    newImgSize = (224,224,0)
    
    if not args.slices is None:
        slices = [int(s) for s in args.slices.split(',')]
    else:
        slices = None

    print('Lambda: %0.6f' %(lambdaBase))

    ## general params
    subjectListCSV = '../input/subjectList_MRU.csv' 

    ## set folders
    FIDsignal = baseFolder + 'FID_folder_fast/subj'+mrn+'_FID_all_avrg_signal.mat'
    FIDRespsignal = baseFolder + 'FID_folder_fast/subj'+mrn+'_FID_respPhaseSignal.mat'
    masksFolder = os.path.dirname(referenceFileSCAN) + '/' 
    
    output_folder = baseFolder + args.output_folder
    print(output_folder)
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    tempFolder =  '/tmp/ltiRecon_subj'+ mrn + '/'
    if not os.path.exists(tempFolder):
        os.makedirs(tempFolder)

    ## set output and temporary files
    coilprofilePath = baseFolder + '/coil_profile_N448/subj'+mrn+'_coilprofile_N448.mat'
    kDataMatFile = tempFolder + '/subj'+mrn+'_kSpaceData.mat'

    ## Load k-space data (use MATLAB script to save to /tmp/subjMRN_kspaceData.mat )
    if (not os.path.exists(kDataMatFile)):
        matFunctionCall = 'addpath(genpath(\'../reconstruction_matlab/ext/\')); try; loadAndSaveKspaceData( {\'%s\'}, \'%s\', %s, \'%s\', \'%s\' );catch; fprintf(\'SOMETHING WENT WRONG WHEN LOADING THE DATA\'); end;exit;' %(patientName,subjectListCSV,str(num_stokes_per_image),kDataMatFile, baseFolder )
        call(['matlab','-nodesktop','-nosplash','-r',matFunctionCall])
        
        
    sliceOverSampling = hdf5storage.loadmat(kDataMatFile, variable_names=['SliceOversampling'])['SliceOversampling']

    ## load respiratory phase
    respData = spio.loadmat(FIDRespsignal)
    
    ## define reconstruction phases. 
    # First corresponds to time, second to respiration.
    # If respiration phase width is 0 or 1, it is ignored
    spkBinWidth = args.stkBinWidth.split(',')
    spkBinWidth = [ float(ph) for ph in spkBinWidth ]
    reconPhases = np.linspace(0,1,respData['Res_Signal'].size).reshape(respData['Res_Signal'].size,1)
    if (len(spkBinWidth) > 1):
        if (spkBinWidth[1] > 0) & (spkBinWidth[1] < 1):
            reconPhases = np.concatenate(( reconPhases, respData['Res_Signal'].reshape(respData['Res_Signal'].size,1) ), axis=-1)

    ## load FIDs
    if args.bulkMotionDetector:
        print( 'Computing outlier spokes' )
        outlierSpokes, outlier_thr, corrFID = detectCorruptedSpokes(FIDsignal, 34, 0 )     # need the number of stokes to create a filtering window. Assuming 34
        spio.savemat( output_folder + 'outlierSpokes.mat', {'outlierSpokes':outlierSpokes,'outlier_thr':outlier_thr,'corrFID':corrFID,'endBaselineIX':0} )
        FID_corr = spio.loadmat(FIDsignal)['FID_corr']
        referenceImgIX = find_reference_img( FID_corr, int(np.round( FID_corr.shape[-1] * spkBinWidth[0] )))   # TODO: this solution might be flawed. It does not account for multiple phases and might have issues with rounding errors
        print(int(np.round( FID_corr.shape[-1] * spkBinWidth[0] )))
        print(referenceImgIX)
    else:
        print( 'Not computing outlier spokes' )
        outlierSpokes = []

    ## set CG optimization and initialize with NUFFT
    print( 'Setting up Conjugate Gradient Solver' )
    cgNUFFT = conjugateGradientRecon( kDataMatFile, coilprofilePath=coilprofilePath, \
                            spkBinWidth=spkBinWidth, spkPhase=reconPhases,\
                            badSpokes=outlierSpokes, slices=slices, num_cores=num_cores )
    imgSize = cgNUFFT.X[0].shape
    T = len(cgNUFFT.X)

    if not args.bulkMotionDetector: # sv407 - set ref image from T if no bulm motion detection 
        referenceImgIX = np.round(T*3/4)
        print(referenceImgIX)

    ## run basic reconstruction (NUFFT) 
    if (args.nufftBool) & ( (not os.path.exists( output_folder +  'nufftRecon_4D.nii.gz')) | (args.force) ):
        print('Computing NUFFT')
        rec_image_NUFFT = cgNUFFT.solve_CG('tikh', lamReg=0.0, alpha_max=10.0, maxIterCG=10, endBaselineIX=0, aortaMask=None)
        graspVol_uint16, flippedGrasp4Dfile, flippedGraspFiles = increaseDimVolume( np.array(rec_image_NUFFT).transpose(1,2,3,0), output_folder +  'nufftRecon', \
                        referenceFileSCAN, flipDim=1, sliceOverSampling=sliceOverSampling, newImgSize=newImgSize, saveIndividualVols=False, rescaleBool=True  )
    else:
        print('Not computing NUFFT')
    
    ## use pyramidal initialization
    if args.pyramInit:
        print('Running pyramidal initialization')
        # pyramParams={'tempRes':[0.2,0.1,0.05,0.01],'lambdas':[lambdaBase,lambdaBase,lambdaBase,lambdaBase]}
        pyramParams={'tempRes':[0.1,0.02],'lambdas':[lambdaBase,lambdaBase]}
        cgNUFFT.pyramidInitialization( pyramParams['tempRes'], pyramParams['lambdas'] )

    ## run basic reconstruction (GRASP) 
    print('Computing GRASP')
    if (not os.path.exists( output_folder +  'graspRecon_4D.nii.gz'))|(args.force):
        for ii in range(2):
            rec_image_GRASP = cgNUFFT.solve_CG('GRASP', lamReg=lambdaBase, alpha_max=10.0, maxIterCG=10, \
                                regMethod=None, endBaselineIX=0, referenceImgIX=referenceImgIX, masksPath=None)
            graspVol_uint16, flippedGrasp4Dfile, flippedGraspFiles = increaseDimVolume( np.array(rec_image_GRASP).transpose(1,2,3,0), output_folder +  'graspRecon_noReg', \
                        referenceFileSCAN, flipDim=1, sliceOverSampling=sliceOverSampling, newImgSize=newImgSize, saveIndividualVols=False, rescaleBool=True  )
        for ii in range(2):
            rec_image_GRASP = cgNUFFT.solve_CG('GRASP', lamReg=lambdaBase, alpha_max=10.0, maxIterCG=10, \
                                    regMethod=args.regMethod, endBaselineIX=0, referenceImgIX=referenceImgIX, masksPath=args.masksPath, sliceOverSampling=sliceOverSampling)
            graspVol_uint16, flippedGrasp4Dfile, flippedGraspFiles = increaseDimVolume( np.array(rec_image_GRASP).transpose(1,2,3,0), output_folder +  'graspRecon', \
                        referenceFileSCAN, flipDim=1, sliceOverSampling=sliceOverSampling, newImgSize=newImgSize, saveIndividualVols=False, rescaleBool=True  )
    else:
        print('\tGRASP already computed')

    ## erase temp folder
    shutil.rmtree(tempFolder, ignore_errors=True)
