# general python imports
import os
import sys
import shutil
from subprocess import call
import argparse
import tempfile
#import matplotlib.pyplot as plt

# math importd
import numpy as np
import nibabel as nib
from scipy.ndimage import morphology as mrp
import scipy.io as spio

import pandas as pd


from joblib import Parallel, delayed
import multiprocessing

sys.path.append('/home/ch199899/Documents/Research/DCE/motion_estimation/')
from sqeuential_registration import parse_affineRegistrationSeries
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
    orderedFiles = orderFiles( goodFiles, [13,7] )

    return orderedFiles

##
#   good luck...
#
def orderFiles( files, posIX ):

    pp = posIX[-1]

    ordering = []
    for ff in files:
        splitSTR = ff.split('_')
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

    num_cores = 20

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
        print maskNames[kidney]


        # if the mask was found, proceed with creating it
        if os.path.exists(masksFolder +  maskNames[kidney]):
            
            dilatedMasks[kidney] = output_folder + maskNames[kidney][:-4]  +'_dilated.nii'
            
            if not os.path.exists(dilatedMasks[kidney]):
                print 'Dilating ' + maskNames[kidney] + ' mask'
                ## load masks
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

            Parallel(n_jobs=num_cores)(delayed(runImageMask)(unmaskedImages[kidney][ii],dilatedMasks[kidney],maskedImages[kidney][ii]) for ii in range( len(unmaskedImages[kidney]) ) )

        else:
            print( masksFolder +  maskNames[kidney] + ' does not exist' )


    ## return path to masked data
    return dilatedMasks, maskedImages, unmaskedImages

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
def register_images_freeForAll( loadFolder, maskedFiles, unmaskedImages, output_folder):


    num_cores = multiprocessing.cpu_count()
    # num_cores = 20

    registrationOutput = range( len(maskedFiles)   )
    
    for ii in range( len(maskedFiles)   ):
        print 'Registering ' + maskedFiles[ii] + ' to everyone'    
        registrationOutput[ii] = Parallel(n_jobs=num_cores)(delayed(rigidRegistration) \
                                    (unmaskedImages[jj], loadFolder + maskedFiles[ii], jj) for jj in range( len(unmaskedImages) ) )


    return registrationOutput

##
#
#
def rigidRegistration( target_Vol, mov_Vol, targetIX ):

    outputVol = output_folder +  os.path.basename(mov_Vol)[:-4] + '_registered_to_vol_' + str(targetIX)

    finalReg = outputVol + '_' + '_rigid'  + '_FinalS.tfm'
    if not os.path.exists( finalReg ):
        ## register volumes
        rigid_reg = outputVol + '_'  + '_rigid'
        print 'Computing: ' + rigid_reg + '_FinalS.tfm'
        if not os.path.exists( rigid_reg + '_FinalS.tfm' ):
            call(['crlBlockMatchingRegistration', '-r', target_Vol, '-f',  mov_Vol, '-o',  rigid_reg \
                                            ,'-p 200', '-t rigid', '--ssi cc', '--inv'])
        finalReg = rigid_reg + '_FinalS.tfm'

        os.remove(rigid_reg + '_FinalS.nrrd')

    return finalReg

##
#
#
def applyRegistration( target_Vol, mov_Vol, regTFM, output_Vol ):
    call([ 'crlResampler2', '-g', target_Vol \
                            , '-i', mov_Vol    \
                            , '-t', regTFM \
                            , '-o', output_Vol ])

##
#   Compute rotation (R) and translation (t) from set of points.
#   Assumes points are in 3D and stacked in columns (N columns for N points)
#
def computeRotationANDTranslation( pin, pout ):

    N = pin.shape[1]    # get number of points

    # get means
    pin_mean = np.mean(pin,axis=1)
    pout_mean = np.mean(pout,axis=1)

    # center points
    X = pin - np.tile( pin_mean, [N,1] ).T
    Y = pout - np.tile( pout_mean, [N,1] ).T

    # compute rotation
    S = X.dot(Y.T)
    U,s,V = np.linalg.svd(S)
    
    R = V.T.dot( np.diag( [1,1,np.linalg.det(V.T.dot(U.T))] ) ).dot(U.T)

    # compute translation
    t = pout_mean - R.dot(pin_mean)

    # move points to new coordinates
    y = R.dot(pin) + np.tile(t,[N,1]).T

    return R, t, y


class reweightedLSQsolver():
        

    ## ---------- MEMBERS ---------- ##
    B = []       # initial rotation of volume i to volume j
    b = []       # initial translation of volume i to volume j
    y = []          # registered position of voxels y[ii] = A[ii]*R + a[ii]
    w = []        # weights of RW-LSQ
    R = []          # position of voxels within mask
    A = []          # rotation result
    a = []          # translation result

    # optimization params
    E0 = 1e-6               # Rw-LSQ regularization
    maxIter = 20            # maximum iterations
    epsilon = 1e-6
    p = 1/2                 # exponent of Lp norm in the objective
    LSQ_epsilon = 1e-6      # epsilon for stopping the LSQ sub-problem
    LSQ_maxIter = 100       # maximum iterations for stopping LSQ sub-problem

    ## ---------- METHODS ---------- ##
    ##
    #   Init function
    #
    def __init__(self, kidneyMaskPath, transformFolder, targetVolIX ):

        # load transform params
        self.load_Transformations( transformFolder )

        # load kidney mask and define position of voxels
        kidneyMask = nib.load( kidneyMaskPath ).get_fdata()

        # compute the position of ROI voxels and store into R
        self.R = np.array( np.where(kidneyMask == 1) )

        # initialize wij
        Nvols = len(self.B)
        self.w = np.ones([Nvols,Nvols])

        # initialize A and a
        for ii in range(Nvols):

            # precompute A
            self.A[ii] = np.eye([3])
            self.a[ii] = np.zeros([3,1])

            # prealocate y
            self.y[ii] = self.B[ii][targetVolIX].dot( self.R ) 

    ##
    #
    #
    #       This is a re-weighted least squares (RW-LSQ) optimization solver for the Lp problem at hand
    #
    #
    def optimize_reweightedLSQ( self ):

        # define
        Nvols = len(self.B)
        Nvox = self.R.shape[1]

        # while no convergence
        kk = 0
        while(True):

            # solve classic LSQ problem
            step_subLSQ = self.solve_subLSQproblem()

            # update Eij and wij
            residual = 0
            for ii in range(Nvols):
                for jj in range(Nvols):
                    # update Eij  >>   Eij = Bij*(Ai*R + ai) + bi -( Aj*R + aj ) (eq. 2)
                    Eij = self.B[ii][jj].dot( self.y[ii] ) + np.tile(self.b[ii][jj],[Nvox,1]).T - self.y[jj]

                    Eij_fro = np.linalg.norm( Eij ,'fro')
                    
                    # update weights >> wij = ( E0 + | Eij |_F )^(p-2)
                    self.w[ii,jj] = ( self.E0 + Eij_fro )**(self.p-2)

                    # update residual
                    residual += Eij_fro / Nvols**2
            
            # restart stepSize
            step_RWLSQ = 0

            # compute rotation and translation matrices and update y
            for ii in range(Nvols):
                y_prev = self.y[ii]
                self.A[ii], self.a[ii], self.y[ii] = computeRotationANDTranslation( self.R, self.y[ii] )

                # update residual
                step_RWLSQ += np.linalg.norm( self.y[ii] - y_prev, 'fro' )

            # report
            print 'Iter: %d. Residual %0.6f. StepSize %0.6f (LSQ step: %0.6f)\n' %( kk, residual, step_RWLSQ, step_subLSQ )

            # check for convergence
            if ( step_RWLSQ < self.epsilon) | ( kk > self.maxIter ):
                print 'Gato Dominguez!'
                return 0

            kk += 1

        return -1

    ##
    #   This method solves the underlying LSQ problem with an alternating directions approach
    #
    def solve_subLSQproblem(self):


        Nvols = len(self.B)

        # while true
        kk = 0
        while(True):

            # residual check for convergence
            residualStep = 0

            # for every y
            for ii in range(Nvols):
                
                # compute big B  TODO: make more efficient using tensor products!
                LHS = np.zeros([3,3])       # left hand side
                RHS = np.zeros([3,1])       # right hand side
                sumIX = [ ix for ix in range(Nvols) if not ix == ii]    # sum for all volumes except for ii
                for jj in sumIX:
                    wBij = self.w[ii,jj] * self.B[ii][jj]
                    RHS += wBij.dot( self.y[jj] + np.tile( self.b[ii][jj], [Nvols,1] ).T )
                    LHS += wBij.T.dot(self.B[ii][jj])

                # solve LSQ subproblem
                y_prev = self.y[ii]
                self.y[ii] = np.linalg.pinv( LHS ).dot( RHS )

                # update residual
                residualStep += np.linalg.norm( self.y[ii] - y_prev, 'fro' )**2

            # check for convergence
            if (residualStep < self.LSQ_epsilon) | ( kk > self.LSQ_maxIter ):
                return  residualStep
            kk += 1

 

    ##
    #   This function takes in a folder with transformations and loads them
    #       It assumes all the transforms in the folder are valid and that there are two levels of ordering (at levels 7 and 13 after splitting the name of the transform)
    #
    #
    def load_Transformations( self, transformFolder ):

        # list all files and filter for .tfm
        tfmFiles = listValidFiles( input_folder, goodTag='.tfm', badTag='', posIX=[7,13] )

        # determine number of time instances (should be sqrt of number of files)
        # and sort transform files into a matrix structure 
        # and load transform and generate rotation matrix + translation vector
        T1 = len( tfmFiles )
        self.B = range(T1)
        self.b = range(T1)
        for t1 in range(T1):
            T2 = len(tfmFiles[t1])
            self.B[t1] = range(T2)
            self.b[t1] = range(T2)
            if not T1 == T2:
                print 'Number of volume tansforms does not match for vol ' + str(tt)

            for t2 in range(T2):
                # sort transform files into a matrix structure 
                AffineTransfrm = parse_affineRegistrationSeries( tfmFiles[t1][t2] )[2]
                # load transform and generate rotation matrix + translation vector
                self.B[t1][t2] = AffineTransfrm[0:3,0:3]
                self.b[t1][t2] = AffineTransfrm[0:3,2]

    def saveAndApplyRegistrationResults(self, saveFolder, mov_Vols, targetIX):

        Nvols = len(self.B)
        
        # for all volumes
        for ii in range(Nvols):
            
            # set name of output
            output_Vol = saveFolder +   os.path.basename(mov_Vols[ii])[:-4] + '_registered_to_vol_' + str(targetIX)

            # save results as TFM transform
            regTFM = output_Vol + '.tfm'
            fo = open(regTFM, 'w+')
            fo.writelines( '#Insight Transform File V1.0\n' )
            fo.writelines( '#Transform 0\n' )
            fo.writelines( 'Parameters:')
            for ll in range(3):
                for jj in range(3):
                    fo.writelines(' %0.8f' %( self.A[ii][jj,ll] ))
            fo.writelines(' %0.8f %0.8f %0.8f' %( self.a[ii][0], self.a[ii][1], self.a[ii][2] ))
            fo.writelines('\n')
            fo.writelines( 'FixedParameters: 0.0000 0.0000 0.0000\n' )
            fo.close()

            # apply TFM transform
            applyRegistration( mov_Vols[targetIX], mov_Vols[ii], regTFM,  output_Vol + '.nrrd' )

    

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

        tempDir = tempfile.mktemp() + '/'
        if not os.path.exists(tempDir):
            os.makedirs(tempDir)


        print 'Loading subject: ' + subj + ' with name: ' + subject_name
        print '\t GRASP folder: ' + input_folder
        print '\t Output folder: ' + output_folder
        print '\t Mask folder: ' + masksFolder


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
            ordering.append( int(splitSTR[7][:-4]) )

        dataFiles = [ goodFiles[ii] for ii in np.argsort(ordering) ]


        ## register all images
        # load FID and determine reference image
        FIDs_all = spio.loadmat(baseFolder + 'FID_folder/s' + subj[1:] + '_FID_all_avrg_signal.mat')
        corrFID = FIDs_all['FID_corr']
        referenceImgIX = find_reference_img(corrFID, num_stokes_per_image)            


        ## create masks
        dilatedMasks, maskedFiles, unmaskedImages = create_movement_masks( masksFolder, input_folder, dataFiles, output_folder + 'tmp/' )

        # for each pair of volumes
        registrationOutput = {'right':[],'left':[]}
        maskNames = {'left':'leftKidneyMask.nii','right':'rightKidneyMask.nii' }
        for kidney in ['left', 'right']:

            ## get masks
            # correct for possible mask names (basically for right kidney... sort of hard coded, I know)
            if not os.path.exists(masksFolder +  maskNames[kidney]):
                if os.path.exists(masksFolder +  maskNames[kidney] + '.gz'):
                    maskNames[kidney] += '.gz'
                elif os.path.exists(masksFolder +  maskNames[kidney][:-4] + 'UpperLower.nii'):
                    maskNames[kidney] = maskNames[kidney][:-4] + 'UpperLower.nii'

            # initialize free-for-all registration
            registrationOutput[kidney] = register_images_freeForAll( input_folder, maskedFiles[kdiney], unmaskedImages, tempDir)

            # refine registration with re-weighted Least Squares Objective
            RWLSQobj = reweightedLSQsolver( maskNames[kidney], tempDir, referenceImgIX )
            RWLSQobj.optimize_reweightedLSQ()
            RWLSQobj.saveAndApplyRegistrationResults(output_folder, maskedFiles['kidney'], referenceImgIX)

    # remove temporary folders
    # os.removedirs(tempDir)