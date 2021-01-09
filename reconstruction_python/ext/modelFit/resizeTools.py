#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  6 09:22:01 2019

@author: ch199899
"""

import os
import numpy as np
import nibabel as nib
from joblib import Parallel, delayed
from scipy.signal import resample

#%%
#       Save the results using the SCANNER format
#
#
#
def increaseDimVolume( volume, outputName, referenceNii=None, flipDim=-1, sliceOverSampling=0, newImgSize=None, saveIndividualVols=False, save4DFile=True, rescaleBool=False, num_cores=20  ):
    
    if referenceNii is None:
        dimNii = newImgSize
        nibReference = None
    else:
        # load reference NII file
        nibReference = nib.load(referenceNii)
        dimNii = nibReference.shape

    # get params
    volDims = volume.shape

    # crop images to the center specified voxels
    if not newImgSize is None:

        Xchange = int((volDims[0] - newImgSize[0] )/2)
        Ychange = int((volDims[1] - newImgSize[1] )/2)

        if (Xchange > 0) & (Ychange > 0):
            yy, xx = np.meshgrid( np.arange( Xchange, volDims[0] - Xchange ), \
                                    np.arange( Ychange, volDims[1] - Ychange )  )
            volume = volume[xx, yy ,:,:]

        else: # extend the image with 0's
            padding = ( (-Xchange, -Xchange), \
                        (-Ychange, -Ychange), \
                        (0,0),(0,0) )
            volume = np.pad( volume, padding, 'constant',constant_values=0  )

    volDims = volume.shape

    # flip desired dimension with FFTSHIFT (if -1 do not do)
    if not flipDim == -1:
        volume = np.fft.fftshift(volume, axes=2)
    
    # evaluate 0 padding needed
    totalslices = np.round( dimNii[2]*(1+sliceOverSampling)) 
    dimDiff = np.round( dimNii[2]*(1+sliceOverSampling) ) - volDims[2] # sila: Consider slice oversampling
    # if np.mod( dimDiff,2 ) == 0:
    #     Zpad1 = int(dimDiff/2)
    #     Zpad2 = Zpad1
    # elif np.mod(dimDiff,2) == 1:
    #     Zpad1 = int(np.floor(dimDiff/2.))
    #     Zpad2 = int(np.abs(dimDiff)-Zpad1) * np.sign(dimDiff)

    # resize Z resolution with FFT
    # if (Zpad1 > 0) & (Zpad2 > 0):
        
    #     # padding = ((0,0),(0,0),(Zpad1,Zpad2),(0,0))
    #     # padded = np.pad( np.fft.fftshift( np.fft.fft(  volume, axis=2) , axes=2 ), padding, 'constant',constant_values=0 )
    #     volume = [ volume[:,:,:,tt] for tt in range(volume.shape[-1]) ]     # dividing problem over time instances to avoid memory issues with long sequences
    #     padding = ((0,0),(0,0),(Zpad1,Zpad2))
    #     def changeDimZ(vol):
    #         padded = np.pad( np.fft.fftshift( np.fft.fft(  vol ,axis=2) , axes=2 ), padding, 'constant',constant_values=0 )
    #         return np.abs( np.fft.ifft( np.fft.ifftshift(  padded ,axes=2) ,axis=2) )
    #     volume = Parallel(n_jobs=num_cores)( delayed(changeDimZ)(volume[tt]) for tt in range(len(volume)) )
    #     volume = np.array(volume).transpose(1,2,3,0)    # rejoining all time instances in single 4D volume
    # elif (Zpad1 < 0) & (Zpad2 < 0):
    #     volume = Parallel(n_jobs=num_cores)( delayed(resample)( volume[:,:,:,tt] , int(totalslices), axis=-1) for tt in range(len(volume)) )
    #     volume = np.array(volume).transpose(1,2,3,0)    # rejoining all time instances in single 4D volume

    orversampledVolume = Parallel(n_jobs=num_cores)( delayed(resample)( volume[:,:,:,tt] , int(totalslices), axis=-1) for tt in range(volume.shape[-1]) )
    orversampledVolume = np.array(orversampledVolume).transpose(1,2,3,0)    # rejoining all time instances in single 4D volume

    volume = orversampledVolume[:,:, \
                        range(int(np.floor((totalslices-dimNii[2])/2.)),int(dimNii[2]+np.floor((totalslices-dimNii[2])/2.))),\
                        :]    #  added by Sila to account for the slice oversampling 


    # keep absolute value of signal only
    volume = np.abs(volume)

    # normalize output to fit in uint16 format
    if rescaleBool:
        maxPixels = np.max(volume.ravel())
        minPixels = np.min(volume.ravel())
        volume_uint16 =  np.array(  np.round( (volume -minPixels) / (maxPixels-minPixels) * 2**15)  , dtype=np.uint16)
	
        print(np.max(volume_uint16))
        print(np.min(volume_uint16))
        # save scaling factors to recover real size
        fo = open(outputName +'_rescaleFactor.txt', 'w+')
        fo.writelines( 'Rescale Factor:\n' )
        fo.writelines( 'maximum value: %0.16f\n' %(maxPixels) )
        fo.writelines( 'minimum value: %0.16f\n' %(minPixels) )
        fo.close()
        
    else:
        volume_uint16 = np.array( np.abs(volume) , dtype=np.float32)

    # save volumes
    indVolumeList = []
    for tt in range(volDims[-1]):

        fileName = outputName + '_' + str(tt) + '.nii.gz'
        indVolumeList.append(fileName)

        if saveIndividualVols:
            if not nibReference is None:
                both_nii = nib.Nifti1Image( volume_uint16[:,:,:,tt], \
                                        nibReference.affine, nibReference.header)
            else:
                both_nii = nib.Nifti1Image( volume_uint16[:,:,:,tt], np.eye(4) )
            both_nii.to_filename( fileName )



    # save 4D volume
    if save4DFile == True:
        fileName = outputName + '_4D.nii.gz'
        if not nibReference is None:
            both_nii = nib.Nifti1Image( volume_uint16, \
                                    nibReference.affine, nibReference.header)
        else:
            both_nii = nib.Nifti1Image( volume_uint16, np.eye(4) )
        both_nii.to_filename( fileName )
    else:
        fileName = None

    return volume_uint16, fileName, indVolumeList
    

       
    
 #%%
#       Save the results using the SCANNER format
#
#
#
def increaseDimMask( mask, outputName, referenceNii=None, flipDim=-1, sliceOverSampling=0, newImgSize=None, saveIndividualVols=False, rescaleBool=False  ):
    
    # get params
    volDims = mask.shape

    # load reference NII file
    nibReference = nib.load(referenceNii)

    # extract sizes
    if not newImgSize is None:
            dimNii = newImgSize
    else:
        dimNii = nibReference.shape
    # print(dimNii)
    ## RESIZE X & Y

    # crop images to the center specified voxels
    if not newImgSize is None:

        Xchange = int((volDims[0] - newImgSize[0] )/2)
        Ychange = int((volDims[1] - newImgSize[1] )/2)

        if (Xchange > 0) & (Ychange > 0):
            yy, xx = np.meshgrid( np.arange( Xchange, volDims[0] - Xchange ), \
                                    np.arange( Ychange, volDims[1] - Ychange )  )
            mask = mask[xx, yy ,:,:]

        else: # extend the image with 0's
            padding = ( (-Xchange, -Xchange), \
                        (-Ychange, -Ychange), \
                        (0,0) )
            
            mask = np.pad( mask, padding, 'constant',constant_values=0  )

    volDims = mask.shape

    ## RESIZE Z

    # consider slice oversampling
    padZ = int(np.floor( volDims[2]*sliceOverSampling/2.) )
    if padZ > 0:     # only apply if really necessary
        padding = ((0,0),(0,0),(padZ,padZ))
        mask = np.pad( mask, padding, 'constant', constant_values=0 )
    
    # change Z dimension
    dimDiff = mask.shape[2] - dimNii[2] # sila: Consider slice oversampling
    if dimDiff > 0:     # only apply if really necessary

        # FFT along Z dim
        Fmask = np.fft.fftshift( np.fft.fft(  mask ,axis=2) , axes=2)

        # truncate Z dim edges
        if np.mod( dimDiff,2 ) == 0:
            Zpad1 = int(dimDiff/2.)
            Zpad2 = Zpad1
        elif np.mod(dimDiff,2) == 1:
            Zpad1 = int(np.floor(dimDiff/2.))
            Zpad2 = int(dimDiff-Zpad1)
        
        Fmask = Fmask[:,:,Zpad1:][:,:,:-Zpad2]
        
        # invert FFT
        mask = np.abs( np.fft.ifft(  np.fft.ifftshift( Fmask ,axes=2) ,axis=2)  )


    ## flip desired dimension with FFTSHIFT
    if not flipDim == -1:
        mask = np.fft.ifftshift(mask, axes=2)
    

    ## normalize output to fit in uint16 format
    if rescaleBool:
        maxPixels = np.max(mask.ravel())
        minPixels = np.min(mask.ravel())
        volume_uint16 =  (mask -minPixels) / (maxPixels-minPixels + 1e-5) 
        volume_uint16[volume_uint16<0.1] = 0
        volume_uint16[volume_uint16>=0.1] = 1
        volume_uint16 =  np.array( np.round( volume_uint16 ) , dtype=np.uint16)
    else:
        volume_uint16 = np.array( mask, dtype=np.float32)


    ## SAVE
    if not os.path.exists(os.path.dirname(outputName)):
        os.makedirs(os.path.dirname(outputName))

    # save volumes
    indVolumeList = []
    if saveIndividualVols:
        for tt in range(volDims[-1]):
            fileName = outputName + '_' + str(tt) + '.nii.gz'
            both_nii = nib.Nifti1Image( volume_uint16[:,:,:,tt], \
                                        nibReference.affine, nibReference.header)
            both_nii.to_filename( fileName )
            indVolumeList.append(fileName)

    # save 4D volume
    fileName = outputName + '.nii.gz'
    both_nii = nib.Nifti1Image( volume_uint16, \
                                nibReference.affine, nibReference.header)
    both_nii.to_filename( fileName )

    return volume_uint16, fileName, indVolumeList
    
