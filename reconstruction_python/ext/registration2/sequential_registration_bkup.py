#!/usr/bin/env python2# -*- coding: utf-8 -*-
"""

@author: jaume Coll-Font <jcollfont@gmail.com>
"""

# general python imports
import os
import sys
import shutil
from subprocess import call
import argparse
from shutil import copyfile
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
#
def registerImages( target_Vol, mov_Vol, mask, outputVol, ix, mode='elastix_dense' ):


    rigid_reg = outputVol + '_'  + '_rigid'
    finalRegImg = rigid_reg + '_FinalS.nii.gz'
    finalRegTFM = rigid_reg + '_FinalS.txt'

    # if not os.path.exists( rigid_reg + '_FinalS.tfm' ):
    if mode == 'elastix_dense':
        elastixFolder = os.path.dirname(outputVol) + '/out%d/' %(ix) 
    
        if not os.path.exists( elastixFolder ):
            os.makedirs(elastixFolder)

        paramFile = '/home/ch199899/Documents/Research/DCE_all/registration/par_pairwise_DCE-ABDOMEN.txt'
        call(['elastix', '-f', target_Vol, '-m', mov_Vol, '-p', paramFile, '-out', elastixFolder ])
        
        copyfile( elastixFolder + '/TransformParameters.0.txt', finalRegTFM)

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
        
    print( 'Done!')

    

    return finalRegTFM


##
#
#
def applyRegistration( target_Vol, mov_Vol, regTFM, output_Vol, ix,  mode='elastix_dense'  ):

    if (mode == 'crl_rigid') | (mode == 'shadab'):
        call([ 'crlResampler2', '-g', target_Vol \
                                , '-i', mov_Vol    \
                                , '-t', regTFM \
                                , '-o', output_Vol ])
    elif  mode == 'elastix_dense' :

        elastixFolder = os.path.dirname(output_Vol) + '/out%d/' %(ix) 
        if not os.path.exists( elastixFolder ):
            os.makedirs(elastixFolder)
        call(['transformix', '-in', mov_Vol, '-out', elastixFolder,'-tp',regTFM])
        if os.path.exists(elastixFolder + '/result.0.nii.gz'):
            copyfile( elastixFolder + '/result.0.nii.gz', output_Vol)
        else:
            copyfile( elastixFolder + '/result.nii.gz', output_Vol)
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
        call([ 'elastix', '-f', imageTarget, '-m', imageTarget, '-t0', transform, '-p',paramFile, '-out', elastixFolder ])

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
#
#
def changeInterpolationMethod( transformParamsFile, outputTransformFile ):

    # eliminate the concatenation of transforms
    fo = open(   transformParamsFile ,'rb')
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

def denseRegistration( target_Vol, mov_Vol, targetIX ):

    outputVol = output_folder +  os.path.basename(mov_Vol)[:-4] + '_registered_to_vol_' + str(targetIX)

    finalReg = outputVol + '_' + '_dense'  + '_FinalS.tfm'
    if not os.path.exists( finalReg ):
        ## register volumes
        rigid_reg = outputVol + '_'  + '_rigid'
        print( 'Computing: ' + rigid_reg + '_FinalS.tfm')
        if not os.path.exists( rigid_reg + '_FinalS.tfm' ):
            call(['crlBlockMatchingRegistration', '-r', target_Vol, '-f',  mov_Vol, '-o',  rigid_reg \
                                            ,'-p 200', '-t rigid', '--ssi cc'])
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
class dciException(Exception):
    pass


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

    print( 'List of Subjects:')
    print( subjects)

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

        print( 'Loading subject: ' + subj + ' with name: ' + subject_name)
        print( '\t GRASP folder: ' + input_folder)
        print( '\t Output folder: ' + output_folder)
        print( '\t Mask folder: ' + masksFolder)

        # create necessary folders
        if not os.path.exists(output_folder):
            os.makedirs(output_folder)
        if not os.path.exists(output_folder + 'tmp/'):
            os.makedirs(output_folder + 'tmp/')

        # detect files in folder
        files = os.listdir( input_folder )
        goodFiles = []
        for ff in files:
            if ff.find(args.tag) > -1:                                  # get files with the 'tag'
                if ( ff.find(args.notag) == -1 ) | ( args.notag == '' ):        # exclude files with the 'notag'
                    goodFiles.append(ff)

        ordering = []
        for ff in goodFiles:
            splitSTR = ff.split('_')
            print( splitSTR)
            ordering.append( int(splitSTR[7][:-4]) )

        dataFiles = [ goodFiles[ii] for ii in np.argsort(ordering) ]

            
        # create dilated masks
        # dilatedMasks, maskedFiles, unmaskedImages = create_movement_masks( masksFolder, input_folder, dataFiles, output_folder + 'tmp/' )

        # # for each pair of volumes
        # registrationOutput = {'right':[],'left':[]}
        # transforms_determinant = {'right':[],'left':[]}
        # transforms_params = {'right':[],'left':[]}
        # residuals_kidneys = {'right':[],'left':[]}
        # image_grad_entropy = {'right':[],'left':[]}
        # image_mean = {'right':[],'left':[]}
        # maskNames = {'left':'leftKidneyMask.nii','right':'rightKidneyMask.nii' }
        # for kidney in ['left', 'right']:

        #     ## get masks
        #     # correct for possible mask names (basically for right kidney... sort of hard coded, I know)
        #     if not os.path.exists(masksFolder +  maskNames[kidney]):
        #         if os.path.exists(masksFolder +  maskNames[kidney] + '.gz'):
        #             maskNames[kidney] += '.gz'
        #         elif os.path.exists(masksFolder +  maskNames[kidney][:-4] + 'UpperLower.nii'):
        #             maskNames[kidney] = maskNames[kidney][:-4] + 'UpperLower.nii'

        ## register all images
        # load FID and determine reference image
        FIDs_all = spio.loadmat(baseFolder + 'FID_folder/s' + subj[1:] + '_FID_all_avrg_signal.mat')
        corrFID = FIDs_all['FID_corr']
        referenceImg = find_reference_img(corrFID, num_stokes_per_image)            

        # apply registration            
        registrationOutput = register_images_toSingle( input_folder, dataFiles, referenceImg, output_folder )


            # ## compute determinant of transforms
            # for ff in registrationOutput[kidney]:
            #     print( ff)
            #     dd, transformParams, A = parse_affineRegistrationSeries( ff )
            #     transforms_determinant[kidney].append(dd)
            #     transforms_params[kidney].append(transformParams)
            #   
            # ## compute image gradient entropy
            # image_grad_entropy[kidney] = []
            # for ff in maskedFiles[kidney]:
            #     img = np.array( nib.load(ff).get_fdata() )
            #     image_grad_entropy[kidney].append( gradient_entropy( [img] ) )

            ## compute residuals and mean (TODO: use tight mask)
        #     residuals_kidneys[kidney], image_mean[kidney] = computeKidneyResiduals( dataFiles[kidney] )

            
        # np.savez( output_folder + 'transforms_determinant.npz', transforms_determinant=transforms_determinant)
        # np.savez( output_folder + 'transforms_params.npz', transforms_params=transforms_params)
        # np.savez( output_folder + 'residuals_kidneys.npz', residuals_kidneys=residuals_kidneys)
        # np.savez( output_folder + 'image_mean_kidneys.npz', image_mean=image_mean)
