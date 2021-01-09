# general python imports
import os
import sys
import shutil
from subprocess import call
import argparse
from shutil import copyfile

# math imports
import numpy as np
import nibabel as nib
import nrrd
import hdf5storage
#import matplotlib.pyplot as plt

# parallel processing imports
from joblib import Parallel, delayed
import multiprocessing


# local imports
sys.path.append('/home/ch199899/Documents/Research/DCE_all/registration/')
from RLTIegistration import runRLTIregistration
from sequential_registration import computeInverseRegTransform, applyRegistration

sys.path.append('/home/ch199899/Documents/Research/DCE_all/modelFit/')
from resizeTools import increaseDimVolume,increaseDimMask

sys.path.append('/home/ch199899/Documents/Research/DCE_all/image_reconstruction/')
from conjugateGradientRecon import conjugateGradientRecon

def invertMasks( maskFiles, outFolder, sliceOverSampling, newImgSize=None, referenceNii=None, rescaleBool=False, flipDim=-1 ):

    inverseMaskFile = []
    for mm in maskFiles:

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
    return flippedLTIFile