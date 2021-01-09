#!/usr/bin/env python2# -*- coding: utf-8 -*-
"""

@author: jaume Coll-Font <jcollfont@gmail.com>
"""

# general python imports
import os
import sys
import shutil
import subprocess
from subprocess import call
import argparse
from shutil import copyfile
import tempfile
#import matplotlib.pyplot as plt

# math importd
import numpy as np
import nibabel as nib
import nrrd
from scipy.ndimage import morphology as mrp
import scipy.io as spio

import transforms3d

import pandas as pd


from joblib import Parallel, delayed
import multiprocessing

# sys.path.append(  os.path.expanduser("~") +  '/links/DCE/motion_estimation/')
# from image_quality_metrics2 import  gradient_entropy


#%%
#
#       This function sequentially registers each image ina sequence to the next.
#       
#   INPUT:
#       - imageList - <T>list - list containing the path to each image.
#       - outputNames - <T>list - list of root names of the results.
#       - maskPath - string - path to the mask used for registration.
#       - badVolumes - <Nbv>int - list of bad volumes.
#       - registrationMode - string - type of registration applied. (refer to function registerImages)
#
#   OUTPUT:
#       - registeredTFM - list - list containing the path to the registered parameters
#
def sequential_registration( imageList, outputNames,  maskPath, badVolumes=[], registrationMode='elastix_dense', num_cores=20 ):

    # select good volumes
    goodVolumes = np.ones([len(imageList)])
    goodVolumes[badVolumes] = 0
    goodVolumes = np.where(goodVolumes)[0]

    T = len(goodVolumes)

    # run registration
    registeredTFM = Parallel(n_jobs=num_cores)(delayed(registerImages) \
                        (imageList[goodVolumes[tt-1]], imageList[goodVolumes[tt]], maskPath, outputNames[goodVolumes[tt]], goodVolumes[tt], mode=registrationMode) for tt in range(1,T ) )

    return registeredTFM
        

##
#
#       This function registers a series of images onto a reference image from the list.
#
#       INPUT:
#           - loadFolder - str - folder where all files are stored
#           - originalImages - str/list - if it is a string, it is the name of the 4D volume containing the data. Otherwise, it should be a list wih every volume to register
#           - output_folder - str - folder where the results will be stored
#           - referenceIX - int/str - (OPTIONAL) volume number (starts at 0) to which all other volumes are registered. If it is a string='sequential', every volume will be registered to its immediate neighbor in the list.
#           - referenceFiles - str/list - (OPTIONAL) if set superseeds refrenceIX. Path to the volume to which all other volumes will be registered.
#           - num_cores - int - (dafault=20)  number of processors to use in parallerl.
#           - tempFolder - str - (optional) temporary folder where to store intermediate results
#           - initialTransform - list - (optional) paths to the initial transforms to use
#           - iterNum - int - (OPTIONAL) modifies the names of the output results. Used in case this function is called iteratively and results at every iteration need to be saved.
#           - defField - str - (OPTIONAL) basename used for the deformation fields. If not specified the deformation fields will not be computed
#
#       OUTPUT:
#           - registrationOutput - list - path to all estimated registration param files
#           - finalRegisteredSequence - str - path to the 4D volume storing the aligned data.
#
def register_images2reference( loadFolder, originalImages, output_folder, referenceIX=0, referenceFiles=None, num_cores=20, tempFolder=None, initialTransform=None, iterNum=0, defField=None ):

    # create temporary directory
    if tempFolder is None:
        tempDir = tempfile.mktemp() + '/'
        os.makedirs(tempDir)
    else:
        tempDir = tempFolder

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # in case the original images is a string pointing at a 4D image
    if isinstance(originalImages, str):
        call(['/home/ch137122/bin/crlConvert4DToN3D','-i', loadFolder + originalImages, '-o', tempDir + os.path.basename(originalImages) ] )
        loadFolder = tempDir
        originalImages = os.listdir( tempDir )
    
    T = len(originalImages)

    # set initial transform
    if initialTransform is None:
        initialTransform = [ None for tt in range(T) ]

    #%% copy images to local temporary folder
    print( 'Saving images to local folder')
    fullImages = []
    for tt in range( T ):
        vol = os.path.basename(originalImages[tt])
        fullImages.append( loadFolder + vol )

    # set deformation field
    if defField:
        defField = [ output_folder + defField.split('.')[0] + '_' + str(tt).zfill(4) + '.nii.gz' for tt in range(T) ]
    else:
        defField = [ None for tt in range(T) ]

    ## set reference images
    if referenceFiles is None:
        if isinstance( referenceIX, (int,float) ):
            referenceFiles = [ fullImages[referenceIX] for tt in range(T) ]
        elif isinstance(referenceIX, str):
            if referenceIX == 'sequential':
                referenceFiles = [ fullImages[tt+1] for tt in range(T) ]   # register every volume to the next
                referenceFiles[T-1] = fullImages[T-1]                      # last volume points at itself
        elif isinstance(referenceIX, (tuple,list)):
            referenceFiles = [ fullImages[tt] for tt in referenceIX ]
    elif isinstance(referenceFiles, str):
        referenceFiles = [ referenceFiles  for tt in range(T) ]
    # elif isinstance(referenceFiles, (tuple,list)):           # no action needed
    
    # prealocate registration names
    registeredImagesName = []
    registeredImages = []
    registrationTransforms = []
    for tt in range(T):
        registeredImagesName.append( tempDir + os.path.basename(fullImages[tt])[:-4] + '_registered_to_vol_' + os.path.basename(referenceFiles[tt][:-4]) )
        rigid_reg = registeredImagesName[-1] + '__rigid'
        registeredImages.append( rigid_reg + '_FinalS.nii.gz' )
        registrationTransforms.append( rigid_reg + '_FinalS.txt' )

    #%% register images to reference volume
    print( 'Computing Registration to reference volume')
    Parallel(n_jobs=num_cores)(delayed(registerImages) \
            (referenceFiles[tt], fullImages[tt], None, registeredImagesName[tt], tt, mode='elastix_dense', initialTransform=initialTransform[tt]) for tt in range( T ) )

    #%% apply  transform to originam images
    print( 'Applying registration to original data')
    Parallel(n_jobs=num_cores)(delayed(applyRegistration) \
            (referenceFiles[tt], fullImages[tt], registrationTransforms[tt], registeredImages[tt], tt, defField=defField[tt]) for tt in range( T ) )

    # join init sequence
    finalRegisteredSequence = os.path.dirname(output_folder) + '/registered_init.nii.gz'
    joinImageSequence( registeredImages, finalRegisteredSequence )

    registrationOutput = []
    for tfmFile in registrationTransforms:
        registrationOutput.append( os.path.dirname(output_folder) + '/TransformParameters.' + str(iterNum) + '.txt' + os.path.basename(tfmFile) )
        try:
            shutil.copy(tfmFile, registrationOutput[-1] )
        except:
            print('Could not copy TFM file: ' + tfmFile + ' to ' + registrationOutput[-1] )
        
    print( 'Gato Dominguez!\n' )
    return registrationOutput, finalRegisteredSequence   

#%% 
#
#   Register 2D images with elastix
#       TODO: assuming FSL format 
#
#
def register2DImagesFrom3Dvols( inputPath, mask=None, referenceVolume=None, outputVolumePath=None, num_cores=40 ):

    # params
    # parFile =  os.path.expanduser("~") + '/Documents/Research/DCE_all/registration2/par_pairwise_rigid-CARDIAC.txt'
    # parFile =  os.path.expanduser("~") + '/Documents/Research/DCE_all/registration2/par_pairwise_deformable-2D.txt'
    parFile =  os.path.expanduser("~") + '/Documents/Research/DCE_all/registration2/par_pairwise_DCE-ABDOMEN.txt'

    # create temp folder
    tempDir = tempfile.mktemp() + '/'
    os.makedirs(tempDir)

    # load volume
    inputFSL = nib.load(inputPath)      # TODO: assuming FSL format 
    inputVolumes = inputFSL.get_fdata()

    # check if input volume is 3D or 4D and set to 4D
    if len( inputVolumes.shape ) == 3:
        inputVolumes = inputVolumes.reshape(inputVolumes.shape + (1,))

    # if reference volume is not provided, set volume 0 of input volume
    if referenceVolume is None:
        referenceVolume = inputVolumes[:,:,:,0]
    elif isinstance(referenceVolume, str):
        referenceVolume = nib.load(referenceVolume).get_fdata()

    # decompose volume into individual slices and time/gradient
    sliceFiles = list(range(inputVolumes.shape[2]))
    for sl in range(inputVolumes.shape[2]):
        sliceFiles[sl] = list(range(inputVolumes.shape[3]))
        for bb in range(inputVolumes.shape[3]):
            sliceFiles[sl][bb] = tempDir + '/vol_sl%d_b%d.nii.gz'%(sl,bb)
            nib.Nifti1Image( inputVolumes[:,:,sl,bb], inputFSL.affine, inputFSL.header).to_filename( sliceFiles[sl][bb] )

    # decompose reference volume into slices
    refFiles = list(range(inputVolumes.shape[2]))
    for sl in range(referenceVolume.shape[2]):
        refFiles[sl] = tempDir + '/refvol_sl%d.nii.gz'%(sl)
        nib.Nifti1Image( referenceVolume[:,:,sl], inputFSL.affine, inputFSL.header).to_filename( refFiles[sl] )

    # prealocate output volume
    outVol4D = np.zeros(inputVolumes.shape)
    print(inputVolumes.shape)

    # apply registration of all slices to a reference slice (TODO: Now this works on a reference volume, but in the future it should be on a slice)
    outSliceFiles = list(range(inputVolumes.shape[2]))
    for sl in range(inputVolumes.shape[2]):
        outSliceFiles[sl] = list(range(inputVolumes.shape[3]))
        Parallel(n_jobs=num_cores)(delayed(registerImages)(refFiles[sl], sliceFiles[sl][bb], mask=mask, outputVol=tempDir + '/out_sl%d_b%d/transform' %(sl,bb), \
                                                            ix=bb, mode='elastix_dense', paramFile=parFile, initialTransform=None ) \
                                                            for bb in range(inputVolumes.shape[3]) )
        # redo 4D volume
        for bb in range(inputVolumes.shape[3]):
            try:
                outVol4D[:,:,sl,bb] = nib.load( tempDir + '/out_sl%d_b%d/out%d/result.0.nii.gz' %(sl,bb,bb)   ).get_fdata()
            except:
                print('Slice %d, volume %d FAILED' %(sl, bb) )

    # save output volume
    if not outputVolumePath is None:
        nib.Nifti1Image( outVol4D, inputFSL.affine, inputFSL.header ).to_filename(outputVolumePath)

    # return registered volume
    return outVol4D


##
#
#
def registerImages( target_Vol, mov_Vol, mask=None, outputVol='', ix=0, mode='elastix_dense', paramFile=None, initialTransform=None ):


    rigid_reg = outputVol + '_'  + '_rigid'
    finalRegImg = rigid_reg + '_FinalS.nii.gz'
    finalRegTFM = rigid_reg + '_FinalS.txt'

    # if not os.path.exists( rigid_reg + '_FinalS.tfm' ):
    if mode == 'elastix_dense':
        elastixFolder = os.path.dirname(outputVol) + '/out%d/' %(ix) 
    
        if not os.path.exists( elastixFolder ):
            os.makedirs(elastixFolder)

        if paramFile is None:
            # paramFile = '/home/ch199899/links/DCE/external/elastix/DownloadedParams/Par0039/par_pairwise/par_real_data/par_pairwise_DCE-ABDOMEN.txt'
            paramFile = '/home/ch199899/Documents/Research/DCE_all/registration2/par_pairwise_DCE-ABDOMEN.txt'

        elastixCommand = ['elastix','-f', target_Vol, '-m', mov_Vol, '-p', paramFile, '-out', elastixFolder ]

        if initialTransform:
            elastixCommand += [ '-t0', initialTransform ]
        if mask:
            elastixCommand += [ '-fMask', mask, '-mMask', mask ]
        
        # print(' '.join(elastixCommand))
        call(elastixCommand, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT )

        try: 
            copyfile( elastixFolder + '/TransformParameters.0.txt', finalRegTFM)
        except:
            print( elastixFolder + '/TransformParameters.0.txt' + ' NOT FOUND. Registration probably failed.' )
        # os.rmdir( elastixFolder )
        
    elif mode == 'shadab':
        ## SHADAB's CODE 
        call( ['/home/ch191070/code/MyCodes/SliceToVolumeRegistration/build/sliceToVol_2' \
                    ,'--NCCmetric'\
                    ,'-m', mov_Vol\
                    ,'-f', target_Vol\
                    ,'-s', mask\
                    ,'-t', finalRegTFM] )

    elif mode == 'crl_rigid':
        print( 'Computing: ' + finalRegImg)
        call(['crlBlockMatchingRegistration', '-r', target_Vol, '-f',  mov_Vol, '-o',  rigid_reg \
                                    ,'-p 200', '-t rigid', '--ssi ssd', '-n 3'])
        
        os.remove(rigid_reg + '_FinalS.nrrd')
        
    # print( 'Done!')

    

    return finalRegTFM


##
#
#
def applyRegistration( target_Vol, mov_Vol, regTFM, output_Vol, ix,  mode='elastix_dense', defField=None  ):

    if (mode == 'crl_rigid') | (mode == 'shadab'):
        call([ 'crlResampler2', '-g', target_Vol \
                                , '-i', mov_Vol    \
                                , '-t', regTFM \
                                , '-o', output_Vol ])

    elif  mode == 'elastix_dense' :

        # create temporary folders
        elastixFolder = os.path.dirname(output_Vol) + '/out%d/' %(ix) 
        if not os.path.exists( elastixFolder ):
            os.makedirs(elastixFolder)

        # setup system call 
        transformixCall = ['transformix', '-in', mov_Vol, '-out', elastixFolder,'-tp',regTFM]
        if defField:
            transformixCall += [ '-def', 'all' ]

        # call transformix
        call( transformixCall, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT )

        # copy results files
        if os.path.exists(elastixFolder + 'result.0.nii.gz'):
            copyfile( elastixFolder + 'result.0.nii.gz', output_Vol)
        else:
            copyfile( elastixFolder + 'result.nii.gz', output_Vol)

        if defField:
            copyfile( elastixFolder + 'deformationField.nii.gz', defField )

        # erase temporary folders
        shutil.rmtree( elastixFolder )

#%%
#
#
def computeInverseRegTransform( imageTarget, transform, inverseTransform, ix=0, finalInterpoaltion=None ):
    
    if not os.path.exists( inverseTransform ):
        
        # create temporary folder for elastix
        elastixFolder = os.path.dirname(inverseTransform) + '/out%d/' %(ix) 
        if not os.path.exists( elastixFolder ):
            os.makedirs(elastixFolder)
        
        paramFile = '/home/ch199899/Documents/Research/DCE_all/registration/par_pairwise_DCE-ABDOMEN.txt'
        call([ 'elastix', '-f', imageTarget, '-m', imageTarget, '-t0', transform, '-p',paramFile, '-out', elastixFolder ]\
                                            , stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT )

        # eliminate the concatenation of transforms
        fo = open(   elastixFolder + '/TransformParameters.0.txt' ,'rb')
        lines = fo.readlines()
        fo.close()

        newLines = []
        for ll in lines:
            if ll.decode('utf8').find('InitialTransformParametersFileName') > -1:
                ll = '(InitialTransformParametersFileName "NoInitialTransform")\n'
                newLines.append(ll)
            elif not finalInterpoaltion is None:
                if ll.decode('utf8').find('(ResampleInterpolator') > -1:
                    ll = '(ResampleInterpolator "FinalNearestNeighborInterpolator")\n'
                    newLines.append(ll)
            else:          
                newLines.append(ll.decode('utf8'))
            

        # write new transform
        fo = open(  inverseTransform ,'w+')
        fo.writelines(newLines)
        fo.close()

        # eliminate elastix folder
        shutil.rmtree( elastixFolder )

    return inverseTransform

#%%
#
#   This function takes an ELASTIX transformation parameter file and concatenates it with another parameter file 'composedTransformFile'.
#
#   If 'composedTransformFile' is left unspecified, the function clears previous composed files and sets 'NoInitialTransform'
#
#
def insertComposedRegistrationParameterFile( originalTransformFile, newTransformFile, composedTransformFile='NoInitialTransform' ):

    # eliminate the concatenation of transforms
    fo = open(  originalTransformFile ,'rb')
    lines = fo.readlines()
    fo.close()

    newLines = []
    for ll in lines:
        if ll.decode('utf8').find('InitialTransformParametersFileName') > -1:
            ll = '(InitialTransformParametersFileName "%s")\n' %(composedTransformFile)
            newLines.append(ll)
        else:          
            newLines.append(ll.decode('utf8'))
        

    # write new transform
    fo = open(  newTransformFile ,'w+')
    fo.writelines(newLines)
    fo.close()


#%%
#
#
#
def changeInterpolationMethod( transformParamsFile, outputTransformFile ):

    # eliminate the concatenation of transforms
    fo = open( transformParamsFile )
    lines = fo.readlines()
    fo.close()

    newLines = []
    for ll in lines:
        # if ll.find('InitialTransformParametersFileName') > -1:
        #     ll = '(InitialTransformParametersFileName "NoInitialTransform")\n'
    
        if ll.find('(ResampleInterpolator') > -1:
            ll = '(ResampleInterpolator "FinalNearestNeighborInterpolator")\n'

        newLines.append(ll)
        

    # write new transform
    fo = open(  outputTransformFile ,'w+')
    fo.writelines(newLines)
    fo.close()

##
#
#   This function aligns pre-existing masks to according to the alignment of two volumes
#
#
#
#
def registerMasks( maskPath, maskVolume, newVolume, newMaskPath ):


    ## register the volume aligned with the mask 'maskVolume' to the new volumes 'newVolume'
    registrationTFM = registerImages( newVolume, maskVolume, mask=None, outputVol=newMaskPath, ix=0, mode='elastix_dense', initialTransform=None )

    ## change interpolation to 
    changeInterpolationMethod( registrationTFM, registrationTFM )

    ## apply registration to masks
    applyRegistration( newVolume, maskPath, registrationTFM, newMaskPath, ix=0,  mode='elastix_dense', defField=None  )

    return registrationTFM, newMaskPath


##
#       
#   This function loads the left and right kidney masks.
#   Then applies a binary opening operation to increase the size of the masks.
#   Finally it applies the mask to the sequence of images.
#
#      This code assumes:
#           - the masks are called leftKidneyMask.nii and rightKidneyMask.nii
#           - the masks are aligned with the target data
#           - availability of crkit for the mathematical operationss
#
#
def create_movement_masks( masksFolder, imgsFolder, anatImgs, output_folder ):

    # for all masks
    dilatedMasks = {'left':[], 'right':[]}
    maskedImages = {'left':[], 'right':[]}
    maskNames = {'left':'leftKidneyMask.nii','right':'rightKidneyMask.nii' }
    unmaskedImages = {'left':[],'right':[] }
    for kidney in ['left', 'right']:
        
        # correct for possible mask names (basically for right kidney... sort of hard coded, I know)
        if not os.path.exists(masksFolder +  maskNames[kidney]):
            if os.path.exists(masksFolder +  maskNames[kidney] + '.gz'):
                maskNames[kidney] += '.gz'
            elif os.path.exists(masksFolder +  maskNames[kidney][:-4] + 'UpperLower.nii'):
                maskNames[kidney] = maskNames[kidney][:-4] + 'UpperLower.nii'
        print( maskNames[kidney])
        # if the mask was found, proceed with creating it
        if os.path.exists(masksFolder +  maskNames[kidney]):
            
            print( 'Dilating ' + maskNames[kidney] + ' mask')

            # compute dilated masks
            dilatedMasks[kidney] = output_folder + maskNames[kidney][:-4]  +'_dilated.nii'

            ## load masks
            # if not os.path.exists( dilatedMasks[kidney] ):
            nii_file = nib.load( masksFolder +  maskNames[kidney] )
            matData = nii_file.get_fdata()

            ## apply binary opening
            matData_filt = mrp.binary_dilation( matData, np.ones([20,20,20]))

            ## save masks
            both_nii = nib.Nifti1Image( np.array(matData_filt, dtype=np.uint16), \
                                        nii_file.affine, nii_file.header)

            both_nii.to_filename(dilatedMasks[kidney] )


            ## multiply masks with data
            for ff in anatImgs:
                unmaskedImages[kidney].append( imgsFolder + ff )
                maskedImages[kidney].append( output_folder + ff[:-4] + '_' + kidney + '_masked.nii' )
                print( 'Masking ' + maskedImages[kidney][-1])
                # if not os.path.exists(maskedImages[kidney][-1]):
                call(['crlImageAlgebra', imgsFolder + ff, 'multiply', dilatedMasks[kidney], maskedImages[kidney][-1] ])
                    

        else:
            print( masksFolder +  maskNames[kidney] + ' does not exist' )


    ## return path to masked data
    return dilatedMasks, maskedImages, unmaskedImages


##
#
#   Loads the tfm file containing all the info about the affine transform and computes its determinant.
#   The objective of this function is to measure "how much" movement there is between two images.
#
#
#
#
def parse_affineRegistrationSeries( tfmFile ):

    fo = open(tfmFile)
    lines = fo.readlines()
    fo.close()

    As = np.eye(4,4)
    eulerParams = []
    nextParams = 0
    for ll in lines:

        # get scaling/rotaion matrices
        if ll.find('AffineTransform_double_3_3') > -1:
            nextParams = 1
        elif ll.find('Euler3DTransform_double_3_3') > -1:
            nextParams = 2
        elif ll.find('VersorRigid3DTransform_double_3_3') > -1:
            nextParams = 2
        elif ll.find('Parameters:') > -1:
            if (nextParams == 1):
                parseLine = ll.split(' ')[1:13]

                A = []
                for nn in parseLine:
                    A.append(float(nn))

                A = np.array(A).reshape(4,3).T
                A = np.concatenate( [A, np.array([[0,0,0,1]]) ],axis=0)
                As = A.dot(As)

                nextParams = 0
            
            elif (nextParams == 2 ):
                parseLine = ll.split(' ')[1:]
                eulerParams = []
                for nn in parseLine:
                    eulerParams.append( float( nn ) )
                
                transformParams = {'translate': eulerParams[3:], 'rotate':eulerParams[:3]}
                nextParams = 0


    sequential_determinant = np.linalg.det(As)

    # transformParams = extractParametersAffineTransform(As)

    return sequential_determinant, transformParams, As

##
#
#       This function saves a new TFM file from a given affine transform
#
#
#
def writeTransformFile(A, tfmFile):

    if A is np.ndarray:
        #%% extract individual parameters
        translate, scale, angles = extractParametersAffineTransform(A)

        #%% prepare lines to write
        lines = []
        lines.append( '#Insight Transform File V1.0\n' )
        lines.append( '#Transform 0\n' )
        lines.append( 'Transform: AffineTransform_double_3_3\n' )
        lines.append( 'Parameters: %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f\n' \
                            %(A[0,0],A[0,1],A[0,2],A[0,3],A[1,0],A[1,1],A[1,2],A[1,3],A[2,0],A[2,1],A[2,2],A[2,3]) )
        lines.append( 'FixedParameters: %0.6f %0.6f %0.6f\n' \
                            %(0,0,0) )
    else:
        lines = []
        lines.append( '#Insight Transform File V1.0\n' )
        lines.append( '#Transform 0\n' )
        lines.append( 'Transform: VersorRigid3DTransform_double_3_3\n' )
        lines.append( 'Parameters: %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f\n' \
                            %(A['rotate'][0],A['rotate'][1],A['rotate'][2], A['translate'][0],A['translate'][1],A['translate'][2]) )
        lines.append( 'FixedParameters: %0.6f %0.6f %0.6f\n' \
                            %(-3.9016716997717253, 33.68466296030961, 0.926670649061478) )

    #%% write results
    fo = open(tfmFile, 'w+')
    fo.writelines(lines)
    fo.close()

    return tfmFile
##
#
#
#
#
def extractParametersAffineTransform(A):

    translate = A[:-1,-1]                       # translation vector

    S = A[:-1,:-1]
    scale = np.linalg.norm(  S ,2,axis=0)       # scaling vector

    R = S / np.tile( scale, [3,1] )             # rotation matrix
    angles = transforms3d.euler.mat2euler(R)    # rotation angles


    return translate, scale, angles

##
#
#
#
#
#
def computeKidneyResiduals( maskedKidneys ):


    # compute residuals
    residuals = [] 
    meanSignal = []
    for ff in range(len(maskedKidneys)-1):
        nibFile = nib.load( maskedKidneys[ff] )
        img1 = nibFile.get_fdata()
        meanSignal.append(  np.mean( img1.ravel()   ) )
        nibFile = nib.load( maskedKidneys[ff+1] )
        img2 = nibFile.get_fdata()
        residuals.append(  np.sum((img1 - img2).ravel()**2) / np.sum(img1.ravel()**2)  )

    if len(maskedKidneys) > 2:
        meanSignal.append(  np.mean( img2.ravel()  ) )
    
    return residuals, meanSignal


##
# determine point with least movement
#
#
# 
def find_reference_img(corrFID,num_stokes_per_image ):

    corrFID_rollSTD = np.zeros( corrFID.shape )
    for zz in range(corrFID.shape[0]):
        series = pd.Series(corrFID[zz,:])
        corrFID_rollSTD[zz,:] = series.rolling(num_stokes_per_image*2).std()

    minimumFID = np.argmin( np.max(corrFID_rollSTD[:,300:],axis=0) ) + 300    # have discarded the first 300 stokes

    referenceImg = int(np.round( minimumFID/ float(num_stokes_per_image) ))

    print( referenceImg)

    return referenceImg
##
#
#
#
#
#
def register_images_toSingle( loadFolder, maskedFiles, targetIX, output_folder):

    target_Vol = loadFolder + maskedFiles[targetIX]

    numCPUs = multiprocessing.cpu_count()

    # registrationOutput = []
    # for ii in range( len(maskedFiles)   ):

    #     # moving volume
    #     mov_Vol = loadFolder + maskedFiles[ii]

    #     # output
    #     outputVol = output_folder +  os.path.basename(mov_Vol)[:-4] + '_registered_to_vol_' + str(targetIX)

    #     print( 'Registering ' + mov_Vol + ' to ' + target_Vol)


    #     # ## register images using Shadab's code
    #     # # if not os.path.exists( outputVol + '.tfm' ):
    #     # call( ['/home/ch191070/code/MyCodes/SliceToVolumeRegistration/build/sliceToVol_2' \
    #     #             ,'--NCCmetric'\
    #     #             ,'-m', mov_Vol\
    #     #             ,'-f', target_Vol\
    #     #             ,'-s', dilatedMasks\
    #     #             ,'-t', outputVol + '.tfm'] )
    #     # registrationOutput.append( outputVol + '.tfm' )

    #     # # if not os.path.exists( outputVol + '.nii' ):
    #     # call([ 'crlResampler2'  , '-g', target_Vol \
    #     #                             , '-i', mov_Vol    \
    #     #                             , '-t', outputVol + '.tfm' \
    #     #                             , '-o', outputVol + '.nii'  ])

    #     ## register volumes DENSE
    #     registrationOutput.append( denseRegistration( target_Vol, mov_Vol ) )
        
    num_cores = 20
    registrationOutput = Parallel(n_jobs=num_cores)(delayed(denseRegistration)(target_Vol, loadFolder + maskedFiles[ii], targetIX) for ii in range( len(maskedFiles) ) )


    ## join in single image
    if (len(registrationOutput) > 0 ):
        joinedFile = output_folder + 'all_registered_to_' + str(targetIX) + '.nrrd'
        # if (not os.path.exists(joinedFile)):
        joinCommand = ['crlConvertN3DTo4D']
        for ff in registrationOutput:
            joinCommand += ['-i', ff[:-4] + '.nrrd']
        joinCommand += ['-o',  joinedFile]
        print( joinCommand)
        call(joinCommand)

    return registrationOutput

##
#
#
#
def denseRegistration( target_Vol, mov_Vol, targetIX ):

    outputVol = output_folder +  os.path.basename(mov_Vol)[:-4] + '_registered_to_vol_' + str(targetIX)

    finalReg = outputVol + '_' + '_dense'  + '_FinalS.tfm'
    if not os.path.exists( finalReg ):
        ## register volumes
        rigid_reg = outputVol + '_'  + '_rigid'
        print( 'Computing: ' + rigid_reg + '_FinalS.tfm')
        if not os.path.exists( rigid_reg + '_FinalS.tfm' ):
            call(['crlBlockMatchingRegistration', '-r', target_Vol, '-f',  mov_Vol, '-o',  rigid_reg \
                                            ,'-p 200', '-t rigid', '--ssi cc'] \
                                            , stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT )
        finalReg = rigid_reg + '_FinalS.tfm'
        
        affine_reg =  outputVol  + '_'  + '_affine'
        print( 'Computing: ' + affine_reg + '_FinalS.tfm')
        if not os.path.exists( affine_reg + '_FinalS.tfm' ):
            call(['crlBlockMatchingRegistration', '-r', target_Vol, '-f',  rigid_reg + '_FinalS.nrrd', '-o',  affine_reg \
                                            , '-i', rigid_reg + '_FinalS.tfm'\
                                            ,'-p 200', '-t affine', '--ssi cc'])
        finalReg = affine_reg + '_FinalS.tfm'
        
        dense_reg = outputVol + '_' + '_dense'
        if not os.path.exists(affine_reg):
            call(['crlBlockMatchingRegistration', '-r', target_Vol, '-f',  affine_reg + '_FinalS.nrrd', '-o',  dense_reg \
                                                , '-i', affine_reg + '_FinalS.tfm'\
                                                ,'-p 200', '-t dense', '--ssi cc'])
        finalReg = dense_reg + '_FinalS.tfm'

        # remove unnecessary volumes
        os.remove(affine_reg + '_FinalS.nrrd')
        os.remove(rigid_reg + '_FinalS.nrrd')

    return finalReg
##
#
#       This function loads the desired files from a given list
#
#
def loadListofFiles( folderPath, filesList, mask ):

    T = len(filesList)

    maskIX = np.where(mask.ravel())[0]

    imgHeaders = range(T)
    Nvx = maskIX.size
    dataVol = np.zeros((Nvx,0))
    for tt in range(T):

        fileFormat = os.path.splitext( filesList[tt] )[-1]
        if fileFormat  == '.nrrd':

            nrrdData = nrrd.read( folderPath + filesList[tt] )
            imgHeaders[tt] = nrrdData[1]
            tmpData = nrrdData[0].transpose(0,2,1)[:,::-1,:].ravel()[maskIX].reshape(maskIX.size,1)
            dataVol = np.concatenate((dataVol, tmpData), axis=1)

        elif (fileFormat  == '.nii') | (fileFormat == '.gz'):

            nii_file = nib.load( folderPath + filesList[tt] )
            tmpData = nii_file.get_fdata().ravel()[maskIX].reshape(maskIX.size,1)
            imgHeaders[tt] = nii_file
            dataVol = np.concatenate((dataVol, tmpData), axis=1)

        else:
            print( 'File format not recognized. Please provide .nii or .nrrd files')
            return dataVol, imgHeaders

        
    return dataVol, imgHeaders

##
#       This function saves a sequence of files provided in a volume to the files in a given lisr
#
def saveVolumeToListOfFiles( dataVol, filesList, imgHeaders, mask ):

    T = len(filesList)

    maskIX = np.where(mask.ravel())[0]
    
    Nx, Ny, Nz = mask.shape

    for tt in range(T):

        individualVolume = dataVol[:,tt].reshape(Nx,Ny,Nz)

        fileFormat = os.path.splitext( filesList[tt] )[-1]
        print( fileFormat)
        if fileFormat  == '.nrrd':
            nrrd.write( filesList[tt], individualVolume[:,::-1,:].transpose(0,2,1), imgHeaders[tt] )

        elif (fileFormat  == '.nii') | (fileFormat == '.gz'):

            both_nii = nib.Nifti1Image( np.array(individualVolume, dtype=np.uint16), \
                                        imgHeaders[tt].affine, imgHeaders[tt].header)

            both_nii.to_filename( filesList[tt] )

        else:
            print( 'File format not recognized. Please provide .nii or .nrrd files')
            return -1

    return 0


##
#
#   This function lists all the files in a folder and selects the files containing the 'goodTag' in their name but no 'badTag'.
#   Then, it sorts the output based on the ennumeration in the names.
#
#   The assumption is that the name is structured with '_' for splitting and that the ordering numbers and in the positions determined by the posIX list. The output will be ordered according to those numbers in the order determined by the list.
#
#
def listValidFiles( input_folder, goodTag='', badTag='', posIX=[] ):

    # detect files in folder
    files = os.listdir( input_folder )
    goodFiles = []
    for ff in files:
        if ff.find(goodTag) | (goodTag=='') > -1:                                  # get files with the 'tag'
            if ( ff.find(badTag) == -1 ) | ( badTag == '' ):        # exclude files with the 'notag'
                goodFiles.append(ff)
    
    # order the files
    orderedFiles = orderFiles( goodFiles, posIX )

    return orderedFiles

##
#   good luck...
#
def orderFiles( files, posIX ):

    pp = posIX[-1]

    ordering = []
    for ff in files:
        splitSTR = ff[:-4].split('_')  # specific to .nii endings..

        if len(splitSTR) >= pp:
            ordering.append( int(splitSTR[pp]) )

    ordering = np.array(ordering)
    uniqueNums = np.unique(ordering)
    orderedFiles = list(range(uniqueNums.size))
    for gg in np.sort(uniqueNums):
        ix = np.where(ordering.ravel()==gg)[0]
        inFiles = [ files[ii] for ii in range(ordering.size) if ordering[ii] == gg ]
        if len(posIX[:-1]) > 0:
            orderedFiles[gg] = orderFiles( inFiles, posIX[:-1] )
        else:
            orderedFiles[gg] = inFiles[0]

    return orderedFiles
    
##
#
#
def joinImageSequence( imageSequence, joinedFile ):
    
    ## join in single image
    if (len(imageSequence) > 0 ):

        # get output folder
        # output_folder = os.path.dirname(imageSequence[0]) + '/'
        
        # create command recursively
        joinCommand = ['/opt/el7/pkgs/crkit/release-current/bin/crlConvertN3DTo4D']
        for ff in imageSequence:
            joinCommand += ['-i', ff]
        joinCommand += ['-o',  joinedFile]
        
        # run
        call(joinCommand)

    return joinedFile       

#%%
#
#   Check if all files in list exist
#
def checkExistListFiles( listFiles ):

    filesExist = True
    for ff in listFiles:
        if not os.path.exists(ff):
            filesExist = False
            break

    return filesExist




#%%
#   
#       generate TFM from rotation and translation params
#       (created to generate synthetic motion for the MEDIA 2020 paper)
#
def generateTFMfromRigidParams( motionParams, newTransforms, refTFM, ResultImagePixelType='unsigned short', imgSize=(224,224,32)):

    for tt in range(len(newTransforms)):

        # create affine matrix
        My = np.array([[ np.cos(motionParams[tt]['rotation'][0]), 0, np.sin(motionParams[tt]['rotation'][0])],\
                        [0, 1, 0],\
                        [-np.sin(motionParams[tt]['rotation'][0]), 0, np.cos(motionParams[tt]['rotation'][0])]])  # rotation around the Y coordinate

        Mz = np.array([[ np.cos(motionParams[tt]['rotation'][1]), -np.sin(motionParams[tt]['rotation'][1]), 0],\
                        [ np.sin(motionParams[tt]['rotation'][1]), np.cos(motionParams[tt]['rotation'][1]), 0],\
                        [ 0, 0, 1] ]) # rotation around the Z coordinate
        
        refAxis = np.array([0,0,1]).reshape(3,1)
        Ux = np.array([[ 0, -refAxis[2,0], -refAxis[1,0]],\
                        [ refAxis[2,0], 0, -refAxis[0,0]],\
                        [ -refAxis[1,0], refAxis[0,0], 0]])
        Ut = np.kron(refAxis,refAxis.T)
        Mr = np.cos(motionParams[tt]['rotation'][2])*np.eye(3) + np.sin(motionParams[tt]['rotation'][2])*Ux + (1-np.cos(motionParams[tt]['rotation'][2] * Ut ) )    # roll around a reference axis

        affineTFM = Mr.dot(Mz.dot(My))
        affineTFM = np.concatenate( (affineTFM, motionParams[tt]['translation'].reshape(3,1)), axis=1 )

        # write new transforms
        fo = open(refTFM[tt])
        origLines = fo.readlines()
        fo.close()

        imageLines = []
        startWrite = False
        center = np.zeros((3,))
        for ll in origLines:

            if (ll.find('Image specific') > -1) | (ll.find('ResampleInterpolator specific') > -1 ) | (ll.find('Resampler specific')>-1):
                startWrite = True
                imageLines.append('\n')
            elif ll.find('// ') > -1:
                startWrite = False

            
            if startWrite:
                
                # modify info from file about image
                if ll.find('(ResultImagePixelType') > -1:
                    ll = '(ResultImagePixelType "%s")\n' %(ResultImagePixelType)
                # if ll.find('(FixedInternalImagePixelType "float")') > -1:
                #     ll = '(FixedInternalImagePixelType "unsigned short")\n'
                # if ll.find('(MovingInternalImagePixelType "float")') > -1:
                #     ll = '(MovingInternalImagePixelType "unsigned short")\n'
                
                if ll.find('(Size') > -1:
                    ll = '(Size'
                    for ii in imgSize:
                        ll += ' %d' %(ii)
                    ll += ')\n'
                if ll.find('FixedImageDimension') > -1:
                    ll = '(FixedImageDimension %d)\n' %(len(imgSize))
                if ll.find('MovingImageDimension') > -1:
                    ll = '(MovingImageDimension %d)\n' %(len(imgSize))
                
                if ll.find('(Spacing ') > -1:
                    temp = ll.split(' ')
                    if len(temp) == len(imgSize):
                        ll = ' '.join(temp)[:-2] + ' 1.0000000000)\n'
                if ll.find('(Origin ') > -1:
                    temp = ll.split(' ')
                    if len(temp) == len(imgSize):
                        ll = ' '.join(temp)[:-2] + ' 0.0000000000)\n'
                if ll.find('(Index ') > -1:
                    temp = ll.split(' ')
                    if len(temp) == len(imgSize):
                        ll = ' '.join(temp)[:-2] + ' 0)\n'
                if ll.find('(Direction ') > -1:
                    temp = ll.split(' ')
                    if len(temp) < len(imgSize)**2:
                        ll = '(Direction '
                        ll += ' '.join(temp[1:4])
                        ll += ' 0.0000000000 '
                        ll += ' '.join(temp[4:7])
                        ll += ' 0.0000000000 '
                        ll += ' '.join(temp[7:10])[:-2]
                        ll += ' 0.0000000000 '
                        ll += ' 0.0000000000 0.0000000000 0.0000000000 1.0000000000)\n'

                imageLines.append(ll)

            if ll.find('(Origin') > -1:
                center = [ float(vv.split(')')[0]) for vv in ll.split(' ')[1:]]
      
        newLines = ['(Transform "AffineTransform")\n']
        if len(imgSize) == 3:
            newLines += ['(NumberOfParameters 12)\n', \
                        '(TransformParameters %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f %0.6f)\n' %( affineTFM[0,0],affineTFM[0,1],affineTFM[0,2],affineTFM[1,0],affineTFM[1,1],affineTFM[1,2],affineTFM[2,0],affineTFM[2,1],affineTFM[2,2],affineTFM[0,3],affineTFM[1,3],affineTFM[2,3] ) ]
        else:
            newLines += ['(NumberOfParameters 20)\n', \
                        '(TransformParameters %0.6f %0.6f %0.6f 0.000000 %0.6f %0.6f %0.6f 0.000000 %0.6f %0.6f %0.6f 0.000000 0.000000 0.000000 0.000000 1.000000 %0.6f %0.6f %0.6f 0.000000)\n' %( affineTFM[0,0],affineTFM[0,1],affineTFM[0,2],affineTFM[1,0],affineTFM[1,1],affineTFM[1,2],affineTFM[2,0],affineTFM[2,1],affineTFM[2,2],affineTFM[0,3],affineTFM[1,3],affineTFM[2,3] ) ]
        newLines += [ '(InitialTransformParametersFileName "NoInitialTransform")\n', \
                    '(UseBinaryFormatForTransformationParameters "false")\n', \
                    '(HowToCombineTransforms "Compose")\n',\
                    '\n// AdvancedAffineTransform specific\n', \
                    '(CenterOfRotationPoint' ]
        for cc in center:
            newLines[-1] += ' %0.6f' %(cc)
        newLines[-1] += ')\n'
        newLines += imageLines

        fo = open( newTransforms[tt], 'w+' )
        fo.writelines(newLines)
        fo.close()

    return newTransforms


#%%
class dciException(Exception):
    pass


## MAIN FUNCTION
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='')
    parser.add_argument('-i', '--input', required=True,
                        help='input to the nrrd/nii 4D volume')
    parser.add_argument('-o', '--output_folder', required=True,
                        help='output folder')
    parser.add_argument('-r', '--refImgIX', required=True,
                        help='Number of reference volume to use as registration (starts at 0). If set to \'sequential\' it will register every volume to the next. If a path to a file is given that data will be used as reference. ') 
    parser.add_argument('-t0', '--initTFM', default=None,
                        help='(OPTIONAL) Initial transform to use in registration') 
    parser.add_argument('-def', '--defField', default=None,
                        help='(OPTIONAL) If given, the deformation field will be computed and saved in a file with this root name' ) 
    parser.add_argument('-m', '--masks', default=None,
                        help='(OPTIONAL) path to the masks to register according to the estimated transforms. Each path should be comma separated. Final interpolation is nearest neighbor' ) 
    parser.add_argument('-p', '--num_cores', default=20,
                        help='Number of cores to use')
    args = parser.parse_args()

    # parse input
    if args.refImgIX.isdigit():
        referenceIX = int(args.refImgIX)
        referenceFiles = None
    elif args.refImgIX == 'sequential':
        referenceIX = args.refImgIX
        referenceFiles = None
    else:
        referenceIX = None
        referenceFiles = args.refImgIX

    num_cores = int(args.num_cores)
    initialTransform = args.initTFM
    defField = args.defField

    # apply registration
    registrationOutput, finalRegisteredSequence = register_images2reference( '', args.input, args.output_folder, referenceIX=referenceIX, referenceFiles=referenceFiles, num_cores=num_cores, tempFolder=None, initialTransform=initialTransform, iterNum=0, defField=defField )

    # # register masks
    # maskPaths = args.m.split(',')
    # for mm in maskPaths:
    #     applyRegistration( mm, mm, registrationOutput[], args.output_folder + os.path.basename( mm )[:-7] + '_regmask', 0,  mode='elastix_dense', defField=None  )
