#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  8 19:02:53 2019

@author: ch199899
"""
import os
import sys

from joblib import Parallel, delayed

# math importd
import numpy as np
import nibabel as nib
import nrrd

sys.path.append('/home/ch199899/Documents/Research/DCE_all/registration/')
from sequential_registration import computeInverseRegTransform, applyRegistration, joinImageSequence, changeInterpolationMethod



def readInMaskAndApplyInverseTransformSeries( maskPath, transformList, outputFolder, referenceIX=0, num_cores=20, inverseImage=None ):

    if inverseImage is None:
        inverseImage = maskPath

    T = len(transformList)

    # prealocate inverseTransformsName
    inverseTransforms = []
    inverseMasks = []
    for tt in range(T):
        inverseTransforms.append(  outputFolder + os.path.basename(transformList[tt])[:-4] + '_inv.txt'  )
        inverseMasks.append( outputFolder + ''.join(os.path.basename(maskPath).split('.')[:-1]) + '_inv_' + str(tt) +  '.nii.gz' )

    # invert transforms
    print( 'Computing inverse transforms')
    Parallel(n_jobs=num_cores)(delayed(computeInverseRegTransform) \
            (maskPath, transformList[tt], inverseTransforms[tt], tt, finalInterpoaltion='FinalNearestNeighborInterpolator') for tt in range( T ) )

    # apply forward transorm to move mask to reference
    referenceMask = outputFolder + ''.join(os.path.basename(maskPath).split('.')[:-1]) + '_reference_'+ str(referenceIX) + '.nii.gz'
    referenceMaskTransform = outputFolder + ''.join(os.path.basename(maskPath).split('.')[:-1]) + '_reference_'+ str(referenceIX) + '_Transform.txt'
    
    print( 'Changing resampler interpolation method to nearest neighbor')
    changeInterpolationMethod( transformList[referenceIX], referenceMaskTransform )
    
    print( 'apply registration to reference mask')
    applyRegistration( maskPath, maskPath, referenceMaskTransform, referenceMask, 0,  mode='elastix_dense'  )

    # apply inverse transforms to the mask at reference position
    print( 'Applying inverse transforms to masks')
    Parallel(n_jobs=num_cores)(delayed(applyRegistration) \
            (referenceMask, referenceMask, inverseTransforms[tt], inverseMasks[tt], tt,  mode='elastix_dense') for tt in range( T ) )

    print( 'Joining Masks')
    joinedInverseMask = outputFolder + ''.join(os.path.basename(maskPath).split('.')[:-1]) + '_inverse.nii.gz'
    joinImageSequence( inverseMasks, joinedInverseMask )
    
    print( 'Loading and returning value')
    return nib.load( joinedInverseMask ).get_fdata()

