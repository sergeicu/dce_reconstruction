# general python imports
import os
import sys
import shutil
import subprocess
from subprocess import call
import argparse
import tempfile
from shutil import copyfile
import glob
#import matplotlib.pyplot as plt

from joblib import Parallel, delayed
import multiprocessing

# math importd
import numpy as np
import nibabel as nib
from scipy.ndimage import morphology as mrp
import scipy.io as spio
import hdf5storage

import pandas as pd

import nrrd

sys.path.append('/home/ch199899/Documents/Research/DCE_all/registration2/')
from sequential_registration import parse_affineRegistrationSeries, writeTransformFile, checkExistListFiles, registerMasks, applyRegistration

sys.path.append('/home/ch199899/Documents/Research/DCE_all/modelFit/')
from rLTI_modelfit import runLTIsysID, runLTISysID4DCE
from resizeTools import increaseDimMask

#%%
#
#
#
#
#
def runRLTIregistration( input_folder, imageFiles, FIDsignal, masksFolder, output_folder, referenceImgIX=None, badVolumes=[], endBaselineIX=0, num_stokes_per_image=34, num_cores=20, tempFolder=None, maxIter=2, referenceNii=None ):

    ## load FID and determine reference image
    if referenceImgIX is None:
        FIDs_all = spio.loadmat(FIDsignal)
        corrFID = FIDs_all['FID_corr']
        referenceImgIX = find_reference_img(corrFID, num_stokes_per_image)            
    
    ## create dilated masks
    kidneyMasks, originalMasks = create_movement_masks( masksFolder,  output_folder, np.ones([20,20,1]), num_cores=num_cores, referenceNii=referenceNii)
    
    ## apply registration  
    registrationOutput, finalGeneratedSequence, finalRegisteredSequence, registrationMaskTFM = \
        register_images_RLTI( input_folder, imageFiles, referenceImgIX, kidneyMasks, output_folder, aifMask=kidneyMasks['aorta'], badVolumes=badVolumes, \
                                endBaselineIX=endBaselineIX, num_cores=num_cores, tempFolder=tempFolder, maxIter=maxIter)

    ## register original masks
    for mm in originalMasks.keys():
        applyRegistration( originalMasks[mm], originalMasks[mm], registrationMaskTFM[mm], output_folder + os.path.basename(originalMasks[mm]), 0,  mode='elastix_dense', defField=None  )

    return registrationOutput, finalGeneratedSequence, finalRegisteredSequence

#%%
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
def create_movement_masks( masksFolder,  output_folder, dilation, num_cores=20, referenceNii=None):

    # prealocate mask path
    dilatedMasks = {'aorta':[]}
    maskNames = {   'left':masksFolder +'leftKidneyMask.nii',\
                    'right':masksFolder +'rightKidneyMask.nii',\
                    'aorta':masksFolder + 'aortaMask.nii' }


    # dilate left and right kidneys adn aorta
    for kidney in ['aorta', 'left', 'right']:
        
        # correct for possible mask names (basically for right kidney... sort of hard coded, I know)
        if not os.path.exists( maskNames[kidney]):
            if os.path.exists(maskNames[kidney] + '.gz'):
                maskNames[kidney] += '.gz'
            elif os.path.exists( maskNames[kidney][:-4] + 'UpperLower.nii'):
                maskNames[kidney] = maskNames[kidney][:-4] + 'UpperLower.nii'
            else:
                filesList = glob.glob( maskNames[kidney][:-4] + '*' )
                if len(filesList) > 0:
                    maskNames[kidney] = filesList[0]

        print( maskNames[kidney])

        # if the mask was found, proceed with creating it
        if os.path.exists( maskNames[kidney]):
            
            dilatedMasks[kidney] = output_folder + os.path.basename(maskNames[kidney]).split('.nii')[0]  +'_dilated.nii.gz'
            
            # if not os.path.exists(dilatedMasks[kidney]):
            print( 'Dilating ' + maskNames[kidney] + ' mask')
            ## load masks
            nii_file = nib.load(  maskNames[kidney] )
            matData = nii_file.get_fdata()

            ## apply binary opening
            matData_filt = mrp.binary_dilation( matData, dilation)
            # matData_filt = mrp.binary_erosion( matData, np.ones([7,7,7]))

            ## save masks 
            if referenceNii is None:
                both_nii = nib.Nifti1Image( np.array(matData_filt, dtype=np.uint16), \
                                            nii_file.affine, nii_file.header)

                both_nii.to_filename( dilatedMasks[kidney] )
            else:   # if reference nii file is provided, save masks according to this size. (adjust so that old mask is cenetered in the middle)
                newMask = increaseDimMask( np.array(matData_filt, dtype=np.uint16), dilatedMasks[kidney][:-7], \
                                            referenceNii=referenceNii, flipDim=-1, sliceOverSampling=0, \
                                            newImgSize=nib.load(referenceNii).shape, saveIndividualVols=False, rescaleBool=True  )


        else:
            print(  maskNames[kidney] + ' does not exist' )

    ## create body and initialization (registration) mask for the rest of the tissue (non aorta or kidneys)
    niiRefFile = nib.load(  dilatedMasks['aorta'] )
    fullMask = np.zeros( niiRefFile.shape, dtype=np.int16)
    for key in dilatedMasks.keys():
        fullMask += np.array( nib.load(  dilatedMasks[key] ).get_fdata(), dtype=np.int16 )

    # create initialization mask
    mskBnd = np.where( fullMask )
    mskBnd = [ np.int16(np.arange( max(0, np.floor(np.min(mskBnd[dd])*0.8)), min(fullMask.shape[dd], np.ceil(np.max(mskBnd[dd])*1.2)) ))  for dd in range(3) ]
    xx,yy = np.meshgrid( mskBnd[0], mskBnd[1] )
    regMask = np.zeros( fullMask.shape, dtype=np.int16)
    regMask[xx,yy,:] = 1
    dilatedMasks['registration'] = output_folder + 'initializationMask.nii.gz'
    tmpNii = nib.Nifti1Image( np.array(regMask, dtype=np.int16), niiRefFile.affine, niiRefFile.header )
    tmpNii.to_filename( dilatedMasks['registration'] )

    # create body mask
    dilatedMasks['body'] = output_folder + 'bodyMask.nii'
    bodyMask = np.ones( niiRefFile.shape, dtype=np.int16)
    bodyMask[ np.where(fullMask > 0) ] = 0
    bodyMask = bodyMask * regMask      # do not consider voxels outside of the init mask TODO: this is being tested now
    #save body mask
    mask_nii = nib.Nifti1Image( np.array(bodyMask, dtype=np.uint16), \
                        niiRefFile.affine, niiRefFile.header)
    mask_nii.to_filename( dilatedMasks['body'] )
          
    print( dilatedMasks )

    ## return path to masked data
    return dilatedMasks, maskNames


def runImageMask(inputImg, mask, outputImg):

    if not os.path.exists(outputImg):
        # print( 'Masking ' + inputImg)
        call(['crlImageAlgebra', inputImg, 'multiply', mask, outputImg ])

    return outputImg

##
#
#
#
#
#
def register_images_RLTI( loadFolder, originalImages, referenceIX, kidneyMasks, output_folder, aifMask=None, badVolumes=[], endBaselineIX=0, num_cores=20, tempFolder=None, maxIter=2):

    # params
    numAtoms = 100
    tau = 1.0
    phiMax = 1e-5
    numClusters = 100
    maxIterLTI = 100
    initParamFile = '/home/ch199899/Documents/Research/DCE_all/registration2/par_pairwise_DCE-ABDOMEN_GRASPrecon448.txt'        # spacing used in B-spline interpolation is 112; much higher than otherwise
    # initParamFile = ' /home/ch199899/Documents/Research/DCE_all/registration2/par_pairwise_DCE-ABDOMEN_coarse.txt'

    # create temporary directory
    if tempFolder is None:
        tempDir = tempfile.mktemp() + '/'
        os.makedirs(tempDir)
    else:
        tempDir = tempFolder
    # tempDir = '/tmp/tmpDxssBK/'   # previousy used folder

    #%% copy images to local temporary folder
    print( 'Saving images to local folder')
    fullImages = []
    if originalImages is None:
        call(['/home/ch137122/bin/crlConvert4DToN3D','-i', loadFolder, '-o', tempDir + 'originalImg.nii.gz'])
        T = nib.load( loadFolder ).shape[-1]
        fullImages = [ tempDir + 'originalImg_'+str(tt).zfill(4)+'.nii.gz' for tt in range(T) ]
        # dataNii = nib.load( loadFolder )
        # data = dataNii.get_fdata()
        # T = dataNii.shape[-1]
        # for tt in range(T):
        #     tmpNii = nib.Nifti1Image( data, dataNii.affine, dataNii.header )
        #     fullImages.append( tempDir + 'originalImg_' + str(tt) + '.nii.gz' )
        #     tmpNii.to_filename( fullImages[-1] )

    else:
        T = len(originalImages)
        for tt in range( T ):
            # original image in vol
            vol = os.path.basename(originalImages[tt])
            fullImages.append( loadFolder + vol )
            # original images registered
            # fullImages.append( tempDir +  os.path.basename(vol) )
            # copyfile( loadFolder + vol, fullImages[-1])     # copy locally

    if aifMask is None:
        aifMask = kidneyMasks['aorta']

    refNii = nib.load(fullImages[0])
    [Nx,Ny,Nz] = refNii.shape

    #%% preload masks
    initMask = kidneyMasks.pop('registration')
    masks = {'aifMask':nib.load(aifMask).get_fdata()}
    bigMask = np.zeros((Nx,Ny,Nz))
    dictKeys = kidneyMasks.keys()
    for kM in dictKeys:
        masks[kM] = nib.load(kidneyMasks[kM]).get_fdata()
        if not kM == 'body':
            bigMask += masks[kM]
    # [Nx,Ny,Nz] = masks['aifMask'].shape
    # mskBnd = np.where( bigMask )
    # mskBnd = [ np.int16(np.arange( max(0, np.floor(np.min(mskBnd[dd])*0.8)), min(bigMask.shape[dd], np.ceil(np.max(mskBnd[dd])*1.2)) ))  for dd in range(3) ]
    # xx,yy = np.meshgrid( mskBnd[0], mskBnd[1] )
    # bigMask[xx,yy,:] = 1
    # initMask = tempDir + 'initializationMask.nii.gz'
    # tmpNii = nib.Nifti1Image( np.array(bigMask, dtype=np.int16), refNii.affine, refNii.header )
    # tmpNii.to_filename( initMask )
          
    ## ----- INITIALIZE ----- ##
    
    # prealocate registration names
    generatedRegisteredImagesName = []
    registeredImagesName = []
    registeredImages = []
    registrationTransforms = []
    generatedRegisteredTransforms = []
    generatedImages = []
    for ff in fullImages:
        registeredImagesName.append( tempDir + os.path.basename(ff)[:-4] + '_registered_to_vol_' + str(referenceIX) )
        
        rigid_reg = registeredImagesName[-1] + '__rigid'
        registeredImages.append( rigid_reg + '_FinalS.nii.gz' )
        registrationTransforms.append( rigid_reg + '_FinalS.txt' )

        generatedImages.append( tempDir + os.path.basename(ff)[:-4] +  '_generated.nii.gz'  )

        generatedRegisteredImagesName.append( tempDir + os.path.basename(ff)[:-4] +  '_generated_registered' )

        generatedRegisteredTransforms.append( generatedRegisteredImagesName[-1] + '__rigid_FinalS.txt' )


    #%% register images to reference volume
    print( 'Computing Initial Registration to reference volume')
    # if not checkExistListFiles( registrationTransforms ):
    Parallel(n_jobs=num_cores)(delayed(rigidRegistration) \
            (fullImages[referenceIX], fullImages[tt], initMask, registeredImagesName[tt], tt, initParamFile ) for tt in range( T ) )
    

    #%% apply  transform to originam images
    print( 'Applying registration to original data')
    # if not checkExistListFiles( registeredImages ):
    Parallel(n_jobs=num_cores)(delayed(applyRegistration) \
            (fullImages[referenceIX], fullImages[tt], registrationTransforms[tt], registeredImages[tt], tt) for tt in range( T ) )

    # join init sequence
    finalGeneratedSequence = None
    finalRegisteredSequence = output_folder + '/' + 'registered_init.nii.gz'
    joinImageSequence( registeredImages, finalRegisteredSequence )

    ## ----- ITERATE ----- ##
    kk = 0
    while( kk < maxIter ):

        ## ----- APPLY rLTI TRUNCATION ----- ##
        print( 'Running LTI model fit')
        # if not checkExistListFiles( generatedImages ):

        print( 'Loading registered volumes')
        imgHeaders = list(range(T))
        inputImg = np.zeros([Nx,Ny,Nz,T])
        for tt in range(T):
            # nrrdData = nrrd.read( registeredImages[tt] )
            # imgHeaders[tt] = nrrdData[1]
            # inputImg[:,:,:,tt] =  nrrdData[0].transpose(0,2,1)[:,::-1,:]
            niiData = nib.load( registeredImages[tt] )
            imgHeaders[tt] = niiData
            inputImg[:,:,:,tt] =  niiData.get_fdata()

        ## get aif (re-compute when working with the aorta and use previously computed otherwise)
        aif = None

        #%% run rLTI
        print( 'Running LTI. This may take a while...')
        ltiImg, aif_rec = runLTISysID4DCE( inputImg, masks, numAtoms=numAtoms, tau=tau, badVolumes=badVolumes, maxIter=maxIterLTI, numClusters=numClusters, nthreads=num_cores, aif=aif, endBaselineIX=endBaselineIX , phiMax=phiMax )

        #%% save results to images
        for tt in range(T):
            # nrrd.write( generatedImages[tt], ltiImg[:,::-1,:,tt].transpose(0,2,1), imgHeaders[tt] )
            both_nii = nib.Nifti1Image( ltiImg[:,:,:,tt], \
                                imgHeaders[tt].affine, imgHeaders[tt].header)
            both_nii.to_filename( generatedImages[tt] )

        aifFile = output_folder + 'arteryInputFunction_iter%d.mat' %(kk)
        spio.savemat(aifFile,{'aif':aif_rec})


        ## ----- REGISTER ORIGINAL IMAGES TO GENERATED IMAGES ----- ##

        #%% register generated (masked) images to full images
        print( 'Registering LTI generated data to full images')
        # if not checkExistListFiles( generatedRegisteredTransforms ):
        Parallel(n_jobs=num_cores)(delayed(rigidRegistration) \
                (generatedImages[tt], fullImages[tt], initMask, generatedRegisteredImagesName[tt], tt ) for tt in range( T ) )

        #%% load transforms
        print( 'Computing inverse transforms')
        stepSizeTrans = 0
        stepSizeRot = 0
        
        ## ----- APPLY INVERSE REGISTERATION TO ORIGINAL IMAGES ----- ##
        registeredImagesName = []
        registrationTransforms = []
        registeredImages = []
        for ff in fullImages:
            registeredImagesName.append( tempDir + os.path.basename(ff)[:-4] + '_registered_to_vol_' + str(referenceIX) )
            rigid_reg = registeredImagesName[-1] + '__rigid' 
            registeredImages.append( rigid_reg + '_FinalS.nii.gz' )
            registrationTransforms.append( rigid_reg + '_FinalS.txt' )

        print( 'Applying registration to original data')
        # if not checkExistListFiles( registeredImages ):
        Parallel(n_jobs=num_cores)(delayed(applyRegistration) \
                (fullImages[referenceIX], fullImages[tt], generatedRegisteredTransforms[tt], registeredImages[tt], tt) for tt in range( T ) )

        ## create 4D images
        print( 'creating 4D volumes for generated and registered images'     )
        # generated
        finalGeneratedSequence = output_folder + '/' + 'generated_iter%d.nii.gz' %(kk)
        joinImageSequence( generatedImages, finalGeneratedSequence )
        # registered
        finalRegisteredSequence = output_folder + '/' + 'registered_iter%d.nii.gz' %(kk)
        joinImageSequence( registeredImages, finalRegisteredSequence )

        # update name of registered images
        generatedRegisteredImagesName = []
        generatedRegisteredTransforms = []
        generatedImages = []
        for ff in fullImages:
            generatedImages.append( tempDir + os.path.basename(ff)[:-4] + '_generated.nii.gz'  )
            generatedRegisteredImagesName.append( tempDir + os.path.basename(ff)[:-4] +  '_generated_registered' )
            generatedRegisteredTransforms.append( generatedRegisteredImagesName[-1] + '__rigid_FinalS.txt')

        ## ----- CHECK CONVERGENCES ----- ##
        print( 'Iter: %d. Step Size: %0.6f, %0.6f' %(kk, stepSizeTrans, stepSizeRot) )
        kk += 1
  
    print( 'Gato Dominguez!\n' )
    
    ## iteration done, set final ouputs and exit
    registrationOutput = []
    for tfmFile in generatedRegisteredTransforms:
        registrationOutput.append( output_folder + '/' + os.path.basename(tfmFile) )
        shutil.copy(tfmFile, registrationOutput[-1] )
        
    ## register masks
    registrationMaskTFM = {}
    for mm in kidneyMasks.keys():
        registrationMaskTFM[mm], newMaskPath = registerMasks( kidneyMasks[mm], fullImages[referenceIX], registeredImages[referenceIX], output_folder + os.path.basename( kidneyMasks[mm] )[:-7] + '_regmask' )

    # delete temporary folder
    shutil.rmtree(tempDir, ignore_errors=True)

    return registrationOutput, finalGeneratedSequence, finalRegisteredSequence, registrationMaskTFM

##
#
#
def rigidRegistration( target_Vol, mov_Vol, mask, outputVol, ix, paramFile=None ):


    rigid_reg = outputVol + '_'  + '_rigid'
    finalRegImg = rigid_reg + '_FinalS.nii.gz'
    finalRegTFM = rigid_reg + '_FinalS.txt'
    # if not os.path.exists( finalRegImg ):
        ## register volumes
    # print( 'Computing: ' + finalRegImg)
        # if not os.path.exists( rigid_reg + '_FinalS.tfm' ):
    # call(['crlBlockMatchingRegistration', '-r', target_Vol, '-f',  mov_Vol, '-o',  rigid_reg \
    #                                 ,'-p 200', '-t rigid', '--ssi ssd', '-n 3'])
        
    # os.remove(rigid_reg + '_FinalS.nrrd')
    # print( 'Done!')
    ## SHADAB's CODE 
    # call( ['/home/ch191070/code/MyCodes/SliceToVolumeRegistration/build/sliceToVol_2' \
    #                 ,'--NCCmetric'\
    #                 ,'-m', mov_Vol\
    #                 ,'-f', target_Vol\
    #                 ,'-s', mask\
    #                 ,'-t', finalRegTFM] )

    elastixFolder = os.path.dirname(outputVol) + '/out%d/' %(ix) 
    if not os.path.exists( elastixFolder ):
        os.makedirs(elastixFolder)

    if paramFile is None:
        # paramFile = '/home/ch199899/links/DCE/external/elastix/DownloadedParams/Par0039/par_pairwise/par_real_data/par_pairwise_DCE-ABDOMEN.txt'
        paramFile = '/home/ch199899/Documents/Research/DCE_all/registration2/par_pairwise_DCE-ABDOMEN.txt'

    elastixCommand = ['elastix', '-f', target_Vol, '-m', mov_Vol, '-p', paramFile, '-out', elastixFolder ]

    if not mask is None:
        elastixCommand += [ '-fMask', mask, '-mMask', mask ]

    call( elastixCommand, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
    
    copyfile( elastixFolder + '/TransformParameters.0.txt', finalRegTFM)
    # os.rmdir( elastixFolder )

    return finalRegImg, finalRegTFM

# ##
# #
# #
# def applyRegistration( target_Vol, mov_Vol, regTFM, output_Vol, ix ):

#     # call([ 'crlResampler2', '-g', target_Vol \
#     #                         , '-i', mov_Vol    \
#     #                         , '-t', regTFM \
#     #                         , '-o', output_Vol ])

#     elastixFolder = os.path.dirname(output_Vol) + '/out%d/' %(ix) 
#     if not os.path.exists( elastixFolder ):
#         os.makedirs(elastixFolder)
#     call(['transformix', '-in', mov_Vol, '-out', elastixFolder,'-tp',regTFM]\
#                                             , stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
#     copyfile( elastixFolder + '/result.nii.gz', output_Vol)

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
        if (ff.find(goodTag) | (goodTag=='') > -1) & ( ff.endswith('.nii.gz') | ff.endswith('.nii') ):                                  # get files with the 'tag'
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

    # retrieve files that have the appropriate tag and do
    ordering = []
    for ff in files:
        if ff.endswith('.nii'):
            splitSTR = ff[:-4].split('_')  # specific to .nii endings..
        elif ff.endswith('.nii.gz'):
            splitSTR = ff[:-7].split('_')  # specific to .nii.gz endings..
        print( splitSTR)
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
    parser.add_argument('-m', '--masks_folder', required=True,
                        help='Masks folder. Must contain a file named aortaMask.nii')
    parser.add_argument('-r', '--refImgIX', required=True,
                        help='Number of reference volume to use as registration (starts at 0).')  
    parser.add_argument('-f', '--corrFIDfile', default=None,
                        help='Path to mat file with the FID metrics for motion detection. (If reference volume is set, this is not necessary)')              
    parser.add_argument('-stk', '--stokes', default='34',
                        help='Number of stokes per image')
    parser.add_argument('-notag', '--notag', default='',
                        help='Non desired tag in the FSL files to load. (e.g. to exclude 4D volume)')
    parser.add_argument('-badvol', '--badvol', default='',
                        help='List of bad volumes separated with a coma')
    parser.add_argument('-strtFrame', '--endBaselineIX', default=0,
                        help='Starting frame of the injection')
    parser.add_argument('-maxIter', '--maxIter', default=2,
                        help='Number of iterations for which the algorithm will be run')                    
    parser.add_argument('-p', '--cores', default=20,
                        help='Number of cores to use')
    args = parser.parse_args()

    num_stokes_per_image = int(args.stokes)
    num_cores = int(args.cores)
    endBaselineIX = int(args.endBaselineIX)
    maxIter = int(args.maxIter)
    referenceImgIX = int(args.refImgIX)

    if not args.badvol == '':
        badVolumes = [int(s) for s in args.badvol.split(',')]
    else:
        badVolumes = []

    # reference folders
    output_folder = args.output_folder
    input_folder = args.input
    masksFolder = args.masks_folder

    # create necessary folders
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    print( 'Loading subject: ' )
    print( '\t GRASP folder: ' + input_folder)
    print( '\t Output folder: ' + output_folder)
    print( '\t Mask folder: ' + masksFolder)

    ## detect valid files in folder and sort them 
    # imageFiles = listValidFiles( input_folder, goodTag=args.tag, badTag=args.notag, posIX=[7] )
    imageFiles = None
    
    ## run RLTI registrationpipeline 
    registrationOutput = runRLTIregistration( input_folder, imageFiles, args.corrFIDfile, masksFolder, output_folder,  \
                                            badVolumes=badVolumes, endBaselineIX=endBaselineIX, num_stokes_per_image=34, referenceImgIX=referenceImgIX, num_cores=num_cores, maxIter=maxIter )

    