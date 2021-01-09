# general python imports
import os
import sys
import shutil
from subprocess import call
import argparse
import tempfile
from shutil import copyfile
#import matplotlib.pyplot as plt

from joblib import Parallel, delayed
import multiprocessing

# math importd
import numpy as np
import nibabel as nib
from scipy.ndimage import morphology as mrp
import scipy.io as spio

import pandas as pd

import nrrd



sys.path.append('/home/ch199899/Documents/Research/DCE/motion_estimation/')
from sequential_registration import parse_affineRegistrationSeries

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

    print referenceImg

    return referenceImg

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
def create_movement_masks( masksFolder, imgsFolder, anatImgs, output_folder , num_cores=20):

    # for all masks
    dilatedMasks = {'left':[], 'right':[]}
    maskedImages = {'left':[], 'right':[]}
    maskNames = {'left':masksFolder +'leftKidneyMask.nii','right':masksFolder +'rightKidneyMask.nii' }
    unmaskedImages = {'left':[],'right':[] }
    for kidney in ['left', 'right']:
        
        # correct for possible mask names (basically for right kidney... sort of hard coded, I know)
        if not os.path.exists( maskNames[kidney]):
            if os.path.exists(maskNames[kidney] + '.gz'):
                maskNames[kidney] += '.gz'
            elif os.path.exists( maskNames[kidney][:-4] + 'UpperLower.nii'):
                maskNames[kidney] = maskNames[kidney][:-4] + 'UpperLower.nii'
        print maskNames[kidney]


        # if the mask was found, proceed with creating it
        if os.path.exists( maskNames[kidney]):
            
            dilatedMasks[kidney] = output_folder + os.path.basename(maskNames[kidney][:-4])  +'_dilated.nii'
            
            if not os.path.exists(dilatedMasks[kidney]):
                print 'Dilating ' + maskNames[kidney] + ' mask'
                ## load masks
                nii_file = nib.load(  maskNames[kidney] )
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

            Parallel(n_jobs=num_cores)(delayed(runImageMask)(unmaskedImages[kidney][ii],dilatedMasks[kidney],maskedImages[kidney][ii]) for ii in range( len(unmaskedImages[kidney]) ) )

        else:
            print(  maskNames[kidney] + ' does not exist' )


    ## return path to masked data
    return dilatedMasks, maskedImages, unmaskedImages, maskNames

def runImageMask(inputImg, mask, outputImg):

    if not os.path.exists(outputImg):
        # print 'Masking ' + inputImg
        call(['crlImageAlgebra', inputImg, 'multiply', mask, outputImg ])

    return outputImg

##
#
#
#
#
#
def register_images_PPCR( loadFolder, maskedImages, referenceImages, referenceIX, kidneyMasks, output_folder, num_cores=20):

    # params
    trunk_SVD = 5
    maxIter = 20
    epsilon = 1e-6

    # create temporary directory
    tempDir = tempfile.mktemp() + '/'
    os.makedirs(tempDir)
    # tempDir = '/tmp/tmpt_sVpQ/'

    # for both kidneys
    registrationOutput = {'left':[], 'right':[]}
    for kidney in ['left','right']:
        
        T = len(unmaskedImages[kidney])
        # initialize data with segmented images
        registeredVolsName = []
        truncatedVolsName = []
        registeredTruncTFMName = []
        for ii in range( T ):
            vol = maskedImages[kidney][ii]
            # original images registered
            registeredVolsName.append( tempDir +  os.path.basename(vol)[:-4] + '_registered_to_vol_' + str(referenceIX) )
            copyfile(vol,registeredVolsName[-1]+'.nii')     # copy locally
            # truncated images
            truncatedVolsName.append( tempDir +  os.path.basename(vol)[:-4] + '_truncated' )
            # truncated & registered imagesoriginal images
            registeredTruncTFMName.append( tempDir +  os.path.basename(vol)[:-4] + '_registeredTruncated_to_vol_' + str(referenceIX) )

        # initialize algorithm with first registration
        print 'Computing Initial Registration'
        registeredVolsName = Parallel(n_jobs=num_cores)(delayed(rigidRegistration) \
                    (referenceImages[kidney][referenceIX], maskedImages[kidney][jj], referenceIX, registeredVolsName[jj]) for jj in range( T ) )
       
        # initialize affine transforms
        AffineTFM = []
        for ff in registeredVolsName:
                tfmFile = ff[0]
                AffineTFM.append( parse_affineRegistrationSeries(tfmFile)[-1] )    # retrieve affine transform from TFM file


        # prealocate kidney indices
        mask = nib.load( kidneyMasks[kidney] ).get_fdata()
        kidneyIX = np.where( mask.ravel() )[0]
        dataMtrx = np.zeros([ kidneyIX.size,T])

        # RUN!
        k = 0
        while(True):
            
            ## ----- PCA TRUNCATION ----- ##

            # read in and join all images into a single sequence
            print 'Loading registered volumes'
            imgHeaders = range(T)
            for tt in range(T):
                nrrdData = nrrd.read( registeredVolsName[tt][0] )
                imgHeaders[tt] = nrrdData[1]
                dataMtrx[:,tt] =  nrrdData[0].ravel()[kidneyIX]
            
            # compute SVD of data
            print 'Computing SVD'
            U,s,V = np.linalg.svd(dataMtrx)

            # project onto truncated subspace
            print 'Projecting Data'
            dataMtrx = U[:,:trunk_SVD].dot( np.diag(s[:trunk_SVD]) ).dot( V[:trunk_SVD,:] )
        
            # save all files
            print 'Saving truncated data'
            for ff in truncatedVolsName:
                tempImg = np.zeros(imgHeaders[ii]['sizes']).ravel()
                tempImg[kidneyIX] = dataMtrx[:,ff]
                nrrd.write( ff, tempImg.reshape(imgHeaders[ii]['sizes']), imgHeaders[ii] )

            ## ----- Registration of Truncated PCA images ----- ##
            print 'Computing registration of truncted data'
            registeredTruncatedTFM = Parallel(n_jobs=num_cores)(delayed(rigidRegistration) \
                                        (referenceImages[kidney][referenceIX], truncatedVolsName[jj], referenceIX, registeredTruncTFMName[jj]) for jj in range( T ) )
            
            ## ----- Apply registration to original images ----- ##
            print 'Applying registration to original data'
            Parallel(n_jobs=num_cores)(delayed(applyRegistration) \
                        (referenceImages[kidney][referenceIX], maskedImages[kidney][jj], registeredTruncatedTFM[jj][1], registeredVolsName[jj][0][:-11]) for jj in range( T ) )

            ## ----- Evaluate convergence ----- ##
            stepSize = 0
            AffineTFM_prev = AffineTFM
            for tt in range(T):
                tfmFile = registeredTruncatedTFM[tt][0]
                AffineTFM[tt] = parse_affineRegistrationSeries(tfmFile)[-1]    # retrieve affine transform from TFM file
                stepSize += np.abs( np.linalg.eig( AffineTFM[tt] - AffineTFM_prev[tt] )[0][0] ) / float(T)

            # check for convergence
            print 'Iter: %d. Step Size: %0.6f' %(k, stepSize) 
            if (k > maxIter) | (stepSize < epsilon):
                print( 'Gato Dominguez!\n' )
                break

            k += 1
        
    
    # delete temporary folder
    # os.removedirs(tempDir)

        
    return registrationOutput

##
#
#
def rigidRegistration( target_Vol, mov_Vol, targetIX, outputVol ):


    rigid_reg = outputVol + '_'  + '_rigid'
    finalRegImg = rigid_reg + '_FinalS.nrrd'
    finalRegTFM = rigid_reg + '_FinalS.tfm'
    # if not os.path.exists( finalRegImg ):
        ## register volumes
    print 'Computing: ' + finalRegImg
        # if not os.path.exists( rigid_reg + '_FinalS.tfm' ):
    call(['crlBlockMatchingRegistration', '-r', target_Vol, '-f',  mov_Vol, '-o',  rigid_reg \
                                    ,'-p 200', '-t rigid', '--ssi cc', '--inv'])
        
        # os.remove(rigid_reg + '_FinalS.nrrd')

    return finalRegImg, finalRegTFM

##
#
#
def applyRegistration( target_Vol, mov_Vol, regTFM, output_Vol ):

    call([ 'crlResampler2', '-g', target_Vol \
                            , '-i', mov_Vol    \
                            , '-t', regTFM \
                            , '-o', output_Vol ])

##
#
#
def joinImageSequence( imageSequence, joinedFile ):
    
    ## join in single image
    if (len(imageSequence) > 0 ):

        # get output folder
        # output_folder = os.path.dirname(imageSequence[0]) + '/'
        
        # create command recursively
        joinCommand = ['crlConvertN3DTo4D']
        for ff in imageSequence:
            joinCommand += ['-i', ff]
        joinCommand += ['-o',  joinedFile]
        
        # run
        call(joinCommand)

    return joinedFile

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
    orderedFiles = orderFiles( goodFiles, [7] )

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
    orderedFiles = range(uniqueNums.size)
    for gg in np.sort(uniqueNums):
        ix = np.where(ordering.ravel()==gg)[0]
        inFiles = [ files[ii] for ii in range(ordering.size) if ordering[ii] == gg ]
        if len(posIX[:-1]) > 0:
            orderedFiles[gg] = orderFiles( inFiles, posIX[:-1] )
        else:
            orderedFiles[gg] = inFiles[0]

    return orderedFiles
        

    

## MAIN FUNCTION
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='')
    parser.add_argument('-i', '--input', required=True,
                        help='input folder to the nrrd sequence')
    parser.add_argument('-tag', '--tag', required=True,
                        help='name tag of the nrrd sequence')
    parser.add_argument('-o', '--output_folder', required=True,
                        help='output folder')
    parser.add_argument('-s', '--subj', default='',
                        help='Name of the subject')
    parser.add_argument('-stk', '--stokes', default='34',
                        help='Number of stokes per image')
    parser.add_argument('-notag', '--notag', default='',
                        help='Non desired tag in the FSL files to load. (e.g. to exclude 4D volume)')
    parser.add_argument('-p', '--cores', default=20,
                        help='Number of cores to use')
    args = parser.parse_args()


    # select subjects
    if args.subj == '':
        files = os.listdir(args.input)
        subjects = []
        for ff in files:
            if ff.find('Subj') > -1:
                subjects.append(ff)
    else:
        subjects = [ args.subj ]

    num_stokes_per_image = int(args.stokes)
    num_cores = int(args.cores)

    print 'List of Subjects:'
    print subjects

    # cycle through all subjects
    for subj in subjects:

        # get subject folder
        baseFolder = args.input
        baseFolder += '/' + subj + '/'
        subject_name = os.listdir(baseFolder)[0] 
        baseFolder += subject_name + '/common_processed/'

        # reference folders
        output_folder = baseFolder + args.output_folder + '/'
        input_folder = baseFolder + '/GRASP_reconstruction_autoBadCoils/'
        masksFolder = baseFolder + '/GRASP_reconstruction_autoBadCoils/'

        # create necessary folders
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        if not os.path.exists(output_folder + 'tmp/'):
            os.makedirs(output_folder + 'tmp/')

        print 'Loading subject: ' + subj + ' with name: ' + subject_name
        print '\t GRASP folder: ' + input_folder
        print '\t Output folder: ' + output_folder
        print '\t Mask folder: ' + masksFolder

        ## detect valid files in folder and sort them 
        dataFiles = listValidFiles( input_folder, goodTag=args.tag, badTag=args.notag, posIX=[7] )

        ## load FID and determine reference image
        FIDs_all = spio.loadmat(baseFolder + 'FID_folder/s' + subj[1:] + '_FID_all_avrg_signal.mat')
        corrFID = FIDs_all['FID_corr']
        referenceImgIX = find_reference_img(corrFID, num_stokes_per_image)            
        
        ## create masks
        dilatedMasks, maskedImages, unmaskedImages, maskNames = create_movement_masks( masksFolder, input_folder, dataFiles, output_folder + 'tmp/', num_cores=num_cores)
        
        ## apply registration            
        registrationOutput = register_images_PPCR( input_folder, maskedImages, unmaskedImages, referenceImgIX, maskNames, output_folder, num_cores=num_cores)