#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on 

@author: Jaume Coll-Font
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
import itertools

# parallel tools
from joblib import Parallel, delayed


# math import
import numpy as np
from scipy.interpolate import interp1d
from scipy.ndimage import morphology as mrp
from pynufft.linalg import nufft_cpu
#from pynufft import NUFFT_cpu

# data import
import nibabel as nib
import nrrd
import hdf5storage
import scipy.io as spio

sys.path.append('/home/ch199899/Documents/Research/DCE_all/registration2/')
from sequential_registration import computeInverseRegTransform, applyRegistration, register_images2reference, insertComposedRegistrationParameterFile
from RLTIegistration import runRLTIregistration

sys.path.append('/home/ch199899/Documents/Research/DCE_all/modelFit/')
from resizeTools import increaseDimVolume,increaseDimMask


##
#
#       This class stores k-space data and implements functions to reconstruct images from it.
#
#
class conjugateGradientRecon():

    ## ----------- MEMBERS ----------- ##

    ## data
    dataPath = ''
    kDataDim = np.zeros((0,))                             # K-space data dimensions
    kData = np.zeros((0,), dtype=np.complex64)            # K-space data
    kPts = np.zeros((0,), dtype=np.complex64)             # position of K-space data
    dcf = np.zeros((0,), dtype=np.complex64)              # density compensation function
    coilprofile = np.zeros((0,), dtype=np.complex64)      # coilprofiles
    X = np.zeros((0,), dtype=np.complex64)                # NUFFT transform of K-space data
    imgDim = np.zeros((0,))                               # recon image dimensions

    ## params
    num_cores = 40      # number of cores to use

    # recon
    numVol = 0          # number of volumes reconstructed
    spkBinWidth = []    # width of each bin used to select spokes in each volume.
    spkPhase = []       # phase (or time) at which each spoke was acquired. 
    spkOfVol = []       # spokes used in volumes
    phaseOfVol = []     # phase of each volume (contains a tuple with the phase in each dim per volume)
    vol2PhaseMap = np.zeros((0,))    # mapping of each volume to the tensor representing the phase maps
    slices = []         # slices to be used in reconstruction
    badSpokes = []      # bad spokes

    # NUFFT
    Nd = (0,0)          # image reconstruction size
    Kd = (0,0)          # oversampling in NUFFT
    Jd = (0,0)          # kernel size in NUFFT

    # CG optimization
    maxIterCG = 10      # number of iterations of Conjugate-Gradient
    gradTollCG = 1e-6   # tolerance of gradient norm in CG. (will be multiplied by the norm of the NUFFT recon)
    eps = 1e-10          # minimum size normalization
    l1smooth = 1e-10     # residual to ensure stability of inverses during the L1 gradient computations

    xNorm = 1.0         # normalization factor for the image residuals/regularizers
    kNorm = 1.0         # normalization factor for the k-space data residuals/regularizers
    
    X = np.zeros((0,), dtype=np.complex64)                 # current reconstructed image
    dX = np.zeros((0,), dtype=np.complex64)                # current step direction
    Fx = np.zeros((0,), dtype=np.complex64)                # current FFT of reconstructed image
    FdX = np.zeros((0,), dtype=np.complex64)               # current FFT of step direction
        
    # line search
    maxlsiter = 150             # maximum number of iterations within the line search
    beta = 0.6                  # scaling factor of the change in step length during line search
    alpha_lin_tol = 0.01        # multiplicative tolerance of the line search
    alpha_max = 1.0             # maximum alpha (initial) in the line search. This will change during computation
    
    # regularization
    ltiModel = np.zeros((0,), dtype=np.float)               # LTI model fit

    # registration  
    regMethod = None                    # type of registration used (if None no registration is applied, speeds computations)
    imageFiles = []                     # list of saved images
    registrationTFM = []                # list of registration transforms
    inverseregistrationTFM = []         # list of inverse registration transforms
    tempDir = ''                        # temporary folder to save results
    regIterNum = 0                      # registration iteration number. Used to keep track of concatenated registrations


    ## ----------- METHODS ----------- ##

    #%%
    #
    #   Class constructor
    #
    #   INPUT:
    #       - dataPath          - str - path to the mat file containing k-space data saved with MATLAB (function: loadAndSaveKspaceData.m )
    #       - coilprofilePath   - str - path to the mat file containing the coil profiles
    #       - slices            - list - list of slices to be reconstructed. If None (default) it uses all available slices.
    #       - spkBinWidth       - double - width of each bin used to select spokes in each volume. If None (default), separates bins into groups of ~34 spokes
    #       - spkPhase          - <1,Spk>double - phase (or time) at which each spoke was acquired. If None (default), it is linear from 0 to 1
    #       - Nd                - <1,2>tuple - X/Y size of the reconstructed images.
    #       - Kd                - <1,2>tuple - X/Y size of the oversampled image in NUFFT.
    #       - Jd                - <1,2>tuple - X/Y size of the kernels used in the NUFFT function.
    #       - badSpokes         - <1,BS>list - list of spokes flagged as outliers
    #       - initX             - list - (OPTIONAL) each element in the list is a 3D volume to initialize the reconstruction.
    #       - num_cores         - int - number of cores to use.
    #       - pyramParams       - dict - (OPTIONAL) if not None sets the pyramidal reconstruction
    #                                   -> 'tempRes': - list - each element in the list corresponds to a spkBinWidth of the pyramidal recon
    #                                   -> 'lambdas': - list - each element in the list corresponds to a lambda to be used in the pyramidal recon
    #                                   Example: pyramParams={'tempRes':[0.2,0.1,0.05,0.01],'lambdas':[10,10,10,10]}
    #
    #
    def __init__(self, dataPath, coilprofilePath=None, slices=None, \
                        spkBinWidth=None, spkPhase=None, \
                        Nd=(448,448), Kd=(int(448*1.5),int(448*1.5)), Jd=(6,6), \
                        badSpokes=[], initX=None, num_cores=20, \
                        pyramParams=None ):

        # set params
        self.Kd = Kd
        self.Nd = Nd
        self.Jd = Jd
        self.num_cores = num_cores

        # load data         (using format inherited from MATLAB load)
        self.dataPath = dataPath
        kspaceData = hdf5storage.loadmat( self.dataPath )
        self.kData = np.complex64( kspaceData['k3n'] )*1e6              ## JCF: scaling by 1^6 to avoid numerical issues
        self.kPts = np.complex64( kspaceData['k_samples'] )
        self.dcf = np.complex64( np.sqrt(  kspaceData['dcf'] ) )       ## JCF: using the squared root of dcf since it will be re-applied multiple times

        # create temporary folder to save results when necessary
        self.tempDir = tempfile.mktemp() + '/'
        print( 'Temp dir: ' + self.tempDir )
        os.makedirs(self.tempDir)


        # select subset of slices for speed (when the slices param is not None)
        if not slices is None:
            self.kData = self.kData[:,:,:,slices]
            self.kPts = self.kPts[:,:,slices]
            self.dcf = self.dcf[:,:,slices]
            self.slices = slices
        else:
            self.slices = np.arange(self.kData.shape[-1])

        self.kDataDim = self.kData.shape    ## # pts x spoke, # channels, # spokes, # slices

        ## PRE-NORMALIZE DATA WITH THE DCF
        # apply dcf normalization and compute square root
        for ch in range(self.kDataDim[1]):
            self.kData[:,ch,:,:] = np.complex64(self.kData[:,ch,:,:] * self.dcf )

        ## PHASE RELATED CALLS AND PREALLOCATIONS
        # set the phase of each spoke if left unspecified
        if spkPhase is None:  
            spkPhase = np.linspace(0,1,self.kDataDim[2]).reshape(self.kDataDim[2],1)
        
        # set the with of each bin in the phase if unspecified
        if spkBinWidth is None:
            spkBinWidth = list(range(spkPhase.shape[1]))
            for dd in range(spkPhase.shape[1]):
                spkBinWidth[dd] = ( np.max(spkPhase[:,dd]) - np.min(spkPhase[:,dd]) ) / 10  # dividing the phase space into 10 parcels per dimension
        
        ## determine spokes in each volume and eliminate bad spokes
        self.defineSpokesOfVol( spkPhase, spkBinWidth, computeDCF=False)
        self.setBadSpokes(badSpokes=badSpokes)  # remove bad spokes from list

        ## LOAD COIL PROFILE
        # load coil profiles
        if coilprofilePath is None:
            self.coilprofile = np.ones(( Nd,Nd, self.kDataDim[1] , self.kDataDim[3] ), dtype=np.complex64)    # if none provided, create flat coilprofile with dims Nd, Nd, NCha, NSli
        else:
            self.coilprofile = np.complex64( hdf5storage.loadmat(coilprofilePath)['coilprofile'] )
            
            self.Nd = self.coilprofile.shape[:2]

            if not slices is None:
                self.coilprofile = self.coilprofile[:,:,:,slices]

        ## INITIAL COMPUTATIONS
        # initialize X by computing the adjoint NUFFT of the recorded data throughtout the entire cycle
        print('Computing initial X')
        if not initX is None:
            self.X = initX
        else:
            self.X = self.iNUFFT()
        
        # self.eps = self.eps*np.max(np.abs( np.array(self.X).ravel() ))
        self.imgDim = self.X[0].shape

        if not pyramParams is None:
            self.pyramidInitialization( pyramParams['tempRes'], pyramParams['lambdas'] )

        # compute normalizations
        self.xNorm = np.prod(self.imgDim)*self.numVol   #np.linalg.norm(np.array(self.X).ravel())**2
        self.kNorm = np.prod(self.kDataDim)             #np.linalg.norm(self.kData.ravel())**2

        # set initial alpha
        self.alpha_max = 1.0#/float(np.prod(self.imgDim)*self.numVol)

        # determine stopping threshold for CG
        # self.gradTollCG = 0
        # for tt in range(self.numVol):
        #     self.gradTollCG += 1e-3*np.linalg.norm( self.X[tt].ravel() ) / self.xNorm
        # self.gradTollCG = 1e-3*np.max(np.abs( np.array(self.X).ravel() ))

        ## REGISTRATION INIT
        # set path to initial registration
        self.registrationTFM = [ os.path.dirname(self.tempDir)  + '/TransformParameters_identity.txt' for tt in range(self.numVol) ]
        self.inverseregistrationTFM = [ os.path.dirname(self.tempDir) + '/TransformParameters_identity.txt' for tt in range(self.numVol) ]

        # create identity transform
        lines = ['(Transform "AffineTransform")\n', \
                    '(NumberOfParameters 12)\n', \
                    '(TransformParameters 1 0 0 0 1 0 0 0 1 0 0 0)\n', \
                    '(InitialTransformParametersFileName "NoInitialTransform")\n', \
                    '(UseBinaryFormatForTransformationParameters "false")\n', \
                    '(HowToCombineTransforms "Compose")\n', \
                    '\n// Image specific\n', \
                    '(FixedImageDimension 3)\n', \
                    '(MovingImageDimension 3)\n', \
                    '(FixedInternalImagePixelType "float")\n', \
                    '(MovingInternalImagePixelType "float")\n', \
                    '(Size %d %d %d)\n' %(self.imgDim[0],self.imgDim[1],self.imgDim[2]), \
                    '(Index 0 0 0)\n', \
                    '(Spacing 1.0 1.0 1.0)\n', \
                    '(Origin 0.0 0.0 0.0)\n', \
                    '(Direction -1.0 0.0 0.0 0.0 -1.0 0.0 0.0 0.0 1.0)\n', \
                    '(UseDirectionCosines "true")\n', \
                    '\n// AdvancedAffineTransform specific\n', \
                    '(CenterOfRotationPoint 0 0 0)\n',
                    '\n// ResampleInterpolator specific\n',
                    '(ResampleInterpolator "FinalBSplineInterpolator")\n',
                    '(FinalBSplineInterpolationOrder 3)\n',
                    '\n// Resampler specific\n',
                    '(Resampler "DefaultResampler")\n',
                    '(DefaultPixelValue 0.000000)\n',
                    '(ResultImageFormat "nii.gz")\n',
                    '(ResultImagePixelType "float")\n']

        fo = open( self.registrationTFM[0], 'w+' )
        fo.writelines(lines)
        fo.close()


    #%%
    #
    #   This function selects the spokes that correspond to each volume according to a parcellation of the phase space (time, respiration, etc.)
    #   Execution of this function results on a new definition of the volumes to be reconstructed and a list of the same length with the phases of each volume.
    #
    #   INPUT:
    #       - phases        - <Spk,D>double - list of phases corresponding to each spoke
    #       - spkBinWidth   - <D,1>double   - width of bins in each phase
    #       - computeDCF    - bool          - if set True, it computes the DCF function from scratch
    #       
    #   OUTPUT/MEMBERS MODIFIED:
    #       - self.numVols       - int   - number of volumes to be reconstructed
    #       - self.spkOfVol      - <numVol,1> list - spokes used to reconstruct each volume
    #       - self.phaseOfVol    - <numVol,D> double - representative phase (mean?) of each volume
    #
    def defineSpokesOfVol( self, phases, spkBinWidth, computeDCF=False ):

        self.spkBinWidth = spkBinWidth
        self.spkPhase = phases

        # set representative dimensions
        D = phases.shape[1]
        Spk = self.kDataDim[2]

        # create phase parcels
        cuttingPoints = list(range(D))
        phasePoints = list(range(D))
        for dd in range(D):
            
            # number of parcels per dim and total number of volumes
            numParcelsDim = int( np.floor( ( np.max(self.spkPhase[:,dd]) - np.min(self.spkPhase[:,dd] )  ) / self.spkBinWidth[dd] ))    # determine number of parcels in phase dimension dd

            # determine cutting points to separate spokes and center phase of each parcel
            cuttingPoints[dd] = [  (pp*self.spkBinWidth[dd] + np.min(self.spkPhase[:,dd]), (pp+1)*self.spkBinWidth[dd] + np.min(self.spkPhase[:,dd])  ) for pp in range(numParcelsDim) ]    # create a tupple with the start and end of each cutting point
            phasePoints[dd] = [  (pp*self.spkBinWidth[dd] + np.min(self.spkPhase[:,dd]) + (pp+1)*self.spkBinWidth[dd] + np.min(self.spkPhase[:,dd]) )/2   for pp in range(numParcelsDim) ]   

        # determine phase of each volume    
        self.phaseOfVol = list(itertools.product(*phasePoints))
        self.numVol = len(self.phaseOfVol)
        volCuttingPoints = list(itertools.product(*cuttingPoints))
        self.vol2PhaseMap = np.arange(self.numVol, dtype=np.int).reshape( [ len(dim) for dim in phasePoints ] )

        # determine spokes per volume
        self.spkOfVol = list(range(self.numVol))
        for tt in range(self.numVol):
            spkIX = np.arange(Spk, dtype=np.int)
            for dd in range(D):
                # get all spokes within range
                spkIX = spkIX[ np.where( (self.spkPhase[spkIX,dd] >= volCuttingPoints[tt][dd][0] ) & (self.spkPhase[spkIX,dd] < volCuttingPoints[tt][dd][1] )  )[0] ]
            self.spkOfVol[tt] = spkIX

        if computeDCF:
            self.computeEmpiricalDCF()

    #%%
    #
    #   This function computes an empirical DCF for every spoke.
    #      TODO: this function will not work properly if there is overlap between the spokes assigned to every volume
    #   
    def computeEmpiricalDCF( self ):

        ## params
        maxiterW = 10000

        # un-apply dcf normalization and compute square root
        for ch in range(self.kDataDim[1]):
            self.kData[:,ch,:,:] = np.complex64( self.kData[:,ch,:,:] / self.dcf )

        def computeSection( timesVec ):

            ## prepare NUFFT object
            NufftObj = NUFFT_cpu()                
            om = self.kPts[:,timesVec,0]
            om  = np.array( [ np.real( om.ravel() )*2*np.pi, np.imag( om.ravel() )*2*np.pi  ] )
            NufftObj.plan( om.T, self.Nd, self.Kd, self.Jd )

            # run forward/inverse operations
            W = self.dcf[:,timesVec,0].ravel()
            for ii in range(maxiterW):
                E = NufftObj.forward( NufftObj.adjoint( W ) )
                W = np.abs(W/E)
            
            return np.sqrt( W ).reshape(self.dcf[:,timesVec,0].shape)     # computing the squared value since it is multiplied twice in the reconstruction

        W = Parallel(n_jobs=self.num_cores)( delayed(computeSection)(self.spkOfVol[tt]) for tt in range(self.numVol) )

        # reassign to DCF variable
        # for tt in range(self.numVol):
        #     self.dcf[:,self.spkOfVol[tt],:] =  np.tile( W[tt].reshape( self.kDataDim[0], len(self.spkOfVol[tt]) ) ,[self.kDataDim[-1], 1,1] ).transpose(1,2,0)

        # average and reassign DCF value (temporary setup)
        print(self.dcf.shape)
        print(np.array(W[0]).shape)
        aveW = np.mean( np.array( [ np.mean(W[tt],axis=1) for tt in range(self.numVol) ]), axis=0 )
        self.dcf = np.tile( aveW, [self.dcf.shape[1],self.dcf.shape[2],1] ).transpose(2,0,1)
        print(self.dcf.shape)

        hdf5storage.savemat('/fileserver/projects6/jcollfont/temp/dcf.mat', {'aveW':aveW})

        # re-apply dcf normalization and compute square root
        for ch in range(self.kDataDim[1]):
            self.kData[:,ch,:,:] = np.complex64(self.kData[:,ch,:,:] * self.dcf )

        return 0


    #%%
    #
    #
    #
    #
    #
    def setBadSpokes(self,badSpokes=[], computeDCF=False):

        self.badSpokes = badSpokes
        # set spokes per volume
        for tt in range(self.numVol):
            self.spkOfVol[tt] = [ spk for spk in self.spkOfVol[tt] if spk not in self.badSpokes ]   # remove bad spokes from list

        # when bad volumes are computed, this must be reapplied
        if computeDCF:
            self.computeEmpiricalDCF()
            

    #%%
    #
    #   Computes inverse NUFFT over all coils and slices. This function parallelizes over time instances.
    #
    #   INPUT:
    #       - kData - <K,Ch,Nspk,Sl>complex - k-space data to be used. If None (default) the function will use the K-space data saved as a class member.
    #
    #   OUTPUT:
    #       - X - <1,numVol>list - list of reconstructed images (over time). Each contains images of size: <Nd[0],Nd[1],Sl>complex
    #
    def iNUFFT( self, kData=None, spkOfVol=None ):

        coilNorm = np.sum( np.abs(self.coilprofile)**2, axis=2 )
        coilDim = self.coilprofile.shape

        if kData is None:
            kData = self.kData
        
        if spkOfVol is None:
            spkOfVol = self.spkOfVol
        numVol = len(spkOfVol)

        def runiNUFFT( timesVec ):

            NufftObj = NUFFT_cpu()
            kDataDim = kData.shape
            
            om = self.kPts[:,timesVec,0]
            om  = [ np.real( om.ravel() )*2*np.pi, np.imag( om.ravel() )*2*np.pi  ] 
            NufftObj.plan(  np.array(om).T, self.Nd, self.Kd, self.Jd )

            xk = list(range(kDataDim[3]))
            for sl in range(kDataDim[3]):      # for all slices
                xk[sl] = np.zeros(self.Nd, dtype=np.complex64)
                if len(timesVec) > 0:
                    for ch in range(kDataDim[1]):      # for all channels 
                        xk[sl] = xk[sl] + np.complex64( NufftObj.adjoint( (kData[:,ch,timesVec,sl] * self.dcf[:,timesVec,sl]).ravel() ) \
                                        * np.conj( self.coilprofile[:,:,ch,sl] ) \
                                        )  / len(timesVec) #/ NufftObj.sn
                    # xk[sl] = xk[sl] * coilDim[0] / np.sqrt( np.prod(coilDim[:2]) ) * np.pi/2. / ( len(timesVec) ) \
                                # / coilNorm[:,:,sl] # normalization over coil profiles  and scaling (JCF: not 100% sold on this one)
                    # xk[sl] = xk[sl] * np.sqrt(np.prod(self.Kd)) 

            # save image in temporary file and apply registration

            return  np.array( xk, dtype=np.complex64).transpose(1,2,0)   # Nd, Nd, Sl
        
        X = Parallel(n_jobs=self.num_cores)( delayed(runiNUFFT)(spkOfVol[tt]) for tt in range(numVol) )
        if not self.regMethod is None:
            X = self.saveAndApplyRegistration( X, self.registrationTFM )
            
        return X

    #%%
    #
    #   Computes forward NUFFT over all coils and slices. This function parallelizes over time instances.
    #
    #   INPUT:
    #       - x - <1, numVol>list - list contatinng the images at each time instane:
    #                                   - <Nd, Nd, Nsli>double - input image to transform to Fourier
    #
    #   OUTPUT:
    #       - kData - <K,Ch,Nspk,Sl>complex - k-space data estimated.
    #
    def NUFFT( self, x ):

        
        def runNUFFT( timesVec, x ):

            NufftObj = NUFFT_cpu()

            imgDim = x.shape

            om = self.kPts[:,timesVec,0]
            om  = [ np.real( om.ravel() )*2*np.pi, np.imag( om.ravel() )*2*np.pi  ] 
            NufftObj.plan(  np.array(om).T, self.Nd, self.Kd, self.Jd ) 

            sk = list(range(self.kDataDim[3]))
            for sl in range(self.kDataDim[3]):      # for all slices
                sk[sl] = list(range(self.kDataDim[1]))
                for ch in range(self.kDataDim[1]):      # for all channels 
                    sk[sl][ch] = np.complex64( NufftObj.forward( self.coilprofile[:,:,ch,sl] * x[:,:,sl] ).reshape( self.kDataDim[0], len(timesVec) )  * self.dcf[:,timesVec,sl] ) #* len(timesVec) #/ np.sqrt( np.prod( imgDim[:2] ) )
            
            return np.array(sk).transpose(2,1,3,0) # K, Ch, Spk, Sl
            
        if not self.regMethod is None:
            x = self.saveAndApplyRegistration( x, self.inverseregistrationTFM )

        return Parallel(n_jobs=self.num_cores)( delayed(runNUFFT)(self.spkOfVol[tt], x[tt] ) for tt in range(self.numVol) )


    #%%
    #
    #   Iteratively solve the regularized LSQ problem to reconstruct the sequence of images X.
    #   This function allows to choose a desired regularization function:
    #       - tikhonov  - 'tikh'  - L2 penalty on the image sequence
    #       - GRASP     - 'GRASP' - L1 penalty on the temporal difference between consecutive images in the sequence
    #       - LTI       - 'LTI'   - L2 penalty on the difference between LTI model fit (saved in ltiModelFile) and the reconstructed image sequence.
    #       
    #   all regularization techniques use lamReg as regularization parameter.
    #
    #   It allows to indicate corrupted spokes (presumably detected with the FID method). Input the spoke index of the corrupted samples as a list in 'badSpokes' when creating the class object
    #
    #
    #       INPUTS:
    #           - regType - str - type of regularization applied
    #           - lamReg - double - amount of regularization applied
    #           - ltiModelFile - str - [required with 'lti' model] path to the file containing the LTI fitted model
    #           - alpha_max - int - [OPTIONAL] starting step length
    #           - X_init - <T>list - [OPTIONAL]  list containing the initialization of recon image. Each entry t contains:
    #                                   - <Nx,Ny,Nz>double - image recon at time instance t
    #
    #
    def solve_CG(self, regType='tikh',  lamReg=0.0, ltiModelFile=None, alpha_max=None, X_init=None, maxIterCG=None, \
                    regMethod=None, masksPath=None, referenceImgIX=None, endBaselineIX=0, sliceOverSampling=0.2  ):   

        ## prepare regularization
        if regType == 'lti':
            self.loadprecomputedLTI(ltiModelFile)

        if not alpha_max is None:
            self.alpha_max = alpha_max

        if not X_init is None:
            self.X = X_init

        if not maxIterCG is None:
            self.maxIterCG = maxIterCG

        if lamReg > 0.0:
            normKdata = [ np.sum(np.abs(self.kData[:,:,self.spkOfVol[tt],:].ravel())**2) for tt in range(self.numVol ) ]
            lamReg = lamReg*np.sum(normKdata)

        ## register
        self.regMethod = regMethod
        if not self.regMethod is None:
            self.registerVolumes( regMethod=self.regMethod, masksPath=masksPath, referenceImgIX=referenceImgIX, endBaselineIX=endBaselineIX, sliceOverSampling=sliceOverSampling)

        ## set initial volume sequence (xk)
        print('Computing Fx')
        self.Fx = self.NUFFT(self.X)

        ## compute initial gradient and step direction (dX)
        print('Precomputing gradient')
        tmpdX = self.computeGradient(regType=regType, lamReg=lamReg )
        dX = [ -1*np.array( tmpdX[tt], dtype=np.complex64) for tt in range(self.numVol) ]
        del tmpdX

        gradKdX_norm = [ np.abs( dX[tt].ravel().conjugate().T.dot(dX[tt].ravel()) ) for tt in range(self.numVol) ]
        gradKdX_norm = np.sum(np.array(gradKdX_norm))
        grad_k_norm = [ ( dX[tt].ravel().conjugate().T.dot(dX[tt].ravel()) ) for tt in range(self.numVol) ]
        grad_k_norm =  np.sum(np.array(grad_k_norm))

        ## while true
        kk = 0
        alpha_prev = 0
        while(True):

            ## line search
            print('Line search')
            alpha, f_min = self.lineSearch( dX, gradKdX_norm, regType=regType, lamReg=lamReg )

            ## update current optimum
            print('Updating xk and computing Fx')
            self.X = [ self.X[tt] + alpha * dX[tt] for tt in range(self.numVol) ]
            
            ## compute step size (plotting purposes)
            stepSize = 0
            for tt in range(self.numVol):
                stepSize += np.linalg.norm( dX[tt].ravel() )

            ## print progress
            print( 'Iter: %d. Cost: %0.6f, step size: %E. Alpha: %E' %(kk, f_min, stepSize, alpha) )
            ## check for convergence
            if ( kk >= self.maxIterCG )| ( alpha*stepSize <  self.gradTollCG ) | (kk >= self.maxIterCG*10):
                print('Gato Dominguez!')

                return self.X

            ## compute new step direction
            print('Computing gradient')
            self.Fx =  self.NUFFT(self.X) 
            grad_k = self.computeGradient( regType=regType, lamReg=lamReg )
            bk = [ grad_k[tt].ravel().conjugate().T.dot(grad_k[tt].ravel()) for tt in range(self.numVol) ]
            bk = np.sum(np.array(bk)) / ( grad_k_norm + self.eps )
            print('bk: %f' %(np.abs(bk)) )
            dX = [ -np.array(  grad_k[tt] + bk * dX[tt], dtype=np.complex64) for tt in range(self.numVol) ]

            grad_k_norm = [ ( grad_k[tt].ravel().conjugate().T.dot(grad_k[tt].ravel()) ) for tt in range(self.numVol) ]
            grad_k_norm =  np.sum(np.array(grad_k_norm))
            gradKdX_norm = [ np.abs( grad_k[tt].ravel().conjugate().T.dot(dX[tt].ravel()) ) for tt in range(self.numVol) ]
            gradKdX_norm = np.sum(np.array(gradKdX_norm))
            del grad_k

            kk +=1
            alpha_prev = alpha

    #%%
    #
    #   This function computes the cost function for the X + alpha*dX.
    #   It assumes that the Fourier transform of the current image sequence X and the step direction dX have been computed!!!
    #
    #   INPUT:
    #       - dX - <1, numVol>list - list contatinng the step direction of the images at each time instane:
    #                                   - <Nd, Nd, Nsli>double - input image to transform to Fourier
    #       - alpha     - double - step length.
    #       - regType   - str - regularization type:
    #                                   - 'tikh' -> tikhonov / l2 norm regularization
    #                                   - 'GRASP' -> GRASP regularization (l1 norm along time)
    #                                   - 'lti' -> l2 difference with LTI model prediction (could be any other prediction)
    #       - lamReg    - double - regularization parameter to be used.
    #
    #   OUTPUT:
    #       - cost - double - optimization function cost.
    #
    def computeCurrentCost(self, dX, alpha, regType='tikh', lamReg=0.0):

        ## main LSQ cost
        LSQcost = Parallel(n_jobs=self.num_cores)( delayed(np.sum)( np.abs(self.Fx[tt] + alpha*self.FdX[tt] - self.kData[:,:,self.spkOfVol[tt],:])**2 ) for tt in range(self.numVol) )
        LSQcost = np.sum(np.array(LSQcost)) 

        ## regularization cost
        if lamReg > 0:
            if regType == 'tikh':
                reguCost = Parallel(n_jobs=self.num_cores)( delayed(np.sum)( np.abs(self.X[tt] + alpha*dX[tt])**2 ) for tt in range(self.numVol) )
                reguCost = np.sum(np.array(reguCost)) 

            elif regType == 'lti':
                reguCost = Parallel(n_jobs=self.num_cores)( delayed(np.sum)( np.abs(self.X[tt] + alpha*dX[tt] - self.ltiModel[:,:,:,tt] )**2 ) for tt in range(self.numVol) )
                reguCost = np.sum(np.array(reguCost)) 

            elif regType == 'GRASP':

                def TVregularization(tt):
                    reguCost = 0
                    ix = np.where(self.vol2PhaseMap == tt )
                    for dd in range( len(ix) ):
                        if ix[dd] > 0:
                            currIX = list(ix)
                            currIX[dd] += -1
                            neighIX = self.vol2PhaseMap[tuple(currIX)][0]
                            reguCost += np.sum(np.sqrt( np.abs( self.X[tt] + alpha*dX[tt] - self.X[neighIX] - alpha*dX[neighIX] )**2  + self.l1smooth ))
                    return lamReg * reguCost

                reguCost = Parallel(n_jobs=self.num_cores)( delayed(TVregularization)( tt ) for tt in range(self.numVol) )
                # reguCost = [ TVregularization( tt ) for tt in range(self.numVol) ]
                reguCost = np.sum(np.array(reguCost)) 

        else:
            reguCost = 0.0

        # print(LSQcost / self.kNorm)
        # print(reguCost/ self.xNorm)
        # print(LSQcost / self.kNorm  +  lamReg * reguCost / self.xNorm)

        return LSQcost / self.kNorm  +  reguCost / self.xNorm

    #%%
    #
    #   This function computes the gradient descent direction.
    #   It is designed to work with various regularization methods, it could be made more general with the introduction of a sub-class for regularization, but might be too cumbersome
    #
    #   INPUT:
    #       - regType   - str - regularization type:
    #                                   - 'tikh' -> tikhonov / l2 norm regularization
    #                                   - 'GRASP' -> GRASP regularization (l1 norm along time)
    #                                   - 'lti' -> l2 difference with LTI model prediction (could be any other prediction)
    #       - lamReg    - double - regularization parameter to be used.
    #
    #   OUTPUT:
    #       - grad - <1,numVol>list - list of gradient images (over time). Each contains images of size: <Nd[0],Nd[1],Sl>complex
    #
    #
    def computeGradient(self, regType='tikh', lamReg=0.0 ):

        ## compute gradient from LSQ
        tempFx = np.zeros( self.kDataDim, dtype=np.complex64 )
        for tt in range(self.numVol):
            tempFx[:,:,self.spkOfVol[tt],:] = 2.*( self.Fx[tt] - self.kData[:,:,self.spkOfVol[tt],:] ) / self.kNorm

        LSQgrad = self.iNUFFT( tempFx )

        ## compute gradient from regularization
        if lamReg > 0:
            if regType == 'tikh':
                for tt in range(self.numVol):
                    LSQgrad[tt]  = LSQgrad[tt] + lamReg* 2.*(self.X[tt]) / self.xNorm
                # reguGrad = [ 2.*( self.X[tt]) for tt in range(self.numVol) ]

            elif regType == 'lti':
                for tt in range(self.numVol):
                    LSQgrad[tt]  = LSQgrad[tt] + lamReg*  2.*( self.X[tt] - self.ltiModel[:,:,:,tt] ) / self.xNorm
                # reguGrad = [ 2.*( self.X[tt] - self.ltiModel[:,:,:,tt] ) for tt in range(self.numVol) ]

            elif regType == 'GRASP':
                
                def TVregularizationGrad( tempGrad, tt):
                    Z1 = np.zeros(self.imgDim, dtype=np.complex64)
                    Z2 = np.zeros(self.imgDim, dtype=np.complex64)
                    ix = np.where(self.vol2PhaseMap == tt )
                    for dd in range( len(ix) ):
                        if ix[dd] > 0:
                            currIX = list(ix)
                            currIX[dd] += -1
                            neighIX = self.vol2PhaseMap[tuple(currIX)][0]
                            Z1 += (self.X[tt] - self.X[neighIX]) / np.sqrt(np.abs(self.X[tt] - self.X[neighIX])**2 + self.l1smooth) 
                        if ix[dd] < self.numVol-1:
                            currIX = list(ix)
                            currIX[dd] += 1
                            neighIX = self.vol2PhaseMap[tuple(currIX)][0]
                            Z2 += (self.X[neighIX] - self.X[tt]) / np.sqrt(np.abs(self.X[neighIX] - self.X[tt])**2 + self.l1smooth) 
                            
                    return tempGrad + lamReg*(Z1 - Z2) / self.xNorm

                LSQgrad = Parallel(n_jobs=self.num_cores)( delayed(TVregularizationGrad)( LSQgrad[tt], tt ) for tt in range(self.numVol) )
                # LSQgrad = [ TVregularizationGrad( tt ) for tt in range(self.numVol) ]

        return LSQgrad
        

    ##
    #
    #   This function loads the LTI fitting results with different storage options (mat or nii files)
    #
    def loadprecomputedLTI(self, ltiModelFile):

        if isinstance( ltiModelFile, list ): 
            T = len(ltiModelFile)
            for tt in range(T):
                nibData = nib.load( ltiModelFile[tt] )
                if tt == 0:
                    self.ltiModel = np.zeros( nibData.get_fdata().shape + (T,))
                self.ltiModel[:,:,:,tt] =  nibData.get_fdata()
        else:
            if ltiModelFile.split('.')[-1] == 'mat':
                self.ltiModel = hdf5storage.loadmat(ltiModelFile)['ltifit']
            elif ltiModelFile.split('.')[-1] == 'nii':
                self.ltiModel = nib.load( ltiModelFile ).get_fdata()

    ##
    #
    #   This function implements a line search algorithm.
    #   It fits a 2nd degree polynomial to the current samples to predict the next alpha.
    #
    def lineSearch2ndOrder(self, dX, gradKdX_norm, regType='tikh', lamReg=0.0):

        # initial costs
        self.FdX = self.NUFFT(dX)
        f = []
        alpha = [0]
        f.append( self.computeCurrentCost( dX, alpha=alpha[0], regType=regType, lamReg=lamReg ) )
        alpha.append( self.alpha_max/2. )
        f.append( self.computeCurrentCost( dX, alpha=alpha[1], regType=regType, lamReg=lamReg ) )
        alpha.append( self.alpha_max )
        f.append( self.computeCurrentCost( dX, alpha=alpha[2], regType=regType, lamReg=lamReg ) )

        ## line search
        lsiter = 0
        f_min = np.min(f)
        while (f_min > f[0] - self.alpha_lin_tol*self.alpha_max*gradKdX_norm) & (lsiter < self.maxlsiter):
            
            ## fit 2nd degree polynomial
            V = np.vander(alpha)
            coefs = np.linalg.pinv( V ).dot( f )

            ## make sure it is not linear
            if coefs[0] > self.eps:
                alpha.append( -coefs[1] / (2.*coefs[0]) )
            else:
                alpha.append( self.alpha_max / self.beta )

            ## evaluate next sample
            f.append( self.computeCurrentCost( dX, alpha=alpha[-1], regType=regType, lamReg=lamReg ) )
            
            ## next iter
            f_min = np.min(f)
            lsiter += 1

        ## save best result
        ix = np.argmin(f)
        alpha_min = alpha[ix]
        print('\tnum iter: %d' %(lsiter))
        
        ## in case of exhausting numer of iterations stop (TODO: probably not necessary)
        if lsiter == self.maxlsiter:
            print('Error - line search ... EXITING')
            self.alpha_max = self.alpha_max * self.beta
            return alpha_min, f_min

        ## update max step length
        if   alpha_min <= self.alpha_max/2.:
            self.alpha_max = alpha_min
        elif alpha_min >= self.alpha_max:
            self.alpha_max = alpha_min / self.beta

        return alpha_min, f_min


    ##
    #
    #   This function implements a line search algorithm.
    #   It fits a 2nd degree polynomial to the current samples to predict the next alpha.
    #
    #   INPUT:
    #       - dX - <1, numVol>list - list contatinng the step direction of the images at each time instane:
    #                                   - <Nd, Nd, Nsli>double - input image to transform to Fourier
    #       - gradKdX_norm  - double - inner product between gradient and step direction.
    #       - regType       - str - regularization type:
    #                                   - 'tikh' -> tikhonov / l2 norm regularization
    #                                   - 'GRASP' -> GRASP regularization (l1 norm along time)
    #                                   - 'lti' -> l2 difference with LTI model prediction (could be any other prediction)
    #       - lamReg    - double - regularization parameter to be used.
    #
    #   OUTPUT:
    #       - alpha_min - double - alpha with minimum cost function found.
    #       - f_min     - double - cost of objective function at minimum lambda
    #
    def lineSearch(self, dX, gradKdX_norm, regType='tikh', lamReg=0.0):

        # initial costs
        self.FdX = self.NUFFT(dX)
        f = []
        alpha = [0]
        f.append( self.computeCurrentCost( dX, alpha=alpha[0], regType=regType, lamReg=lamReg ) )
        alpha.append( self.alpha_max )
        f.append( self.computeCurrentCost( dX, alpha=alpha[1], regType=regType, lamReg=lamReg ) )

        ## line search
        lsiter = 0
        while (f[-1] > f[0] - self.alpha_lin_tol*alpha[-1]*gradKdX_norm) & (lsiter < self.maxlsiter):
            
            ## update alpha
            alpha.append( alpha[-1] * self.beta )

            ## evaluate next sample
            f.append( self.computeCurrentCost( dX, alpha=alpha[-1], regType=regType, lamReg=lamReg ) )
            
            ## next iter
            f_min = np.min(f)
            lsiter += 1

            ## check if algrithm is doing dumb computations
            # print(np.abs( f[-1] - f[0] ))
            if np.abs( f[-1] - f[0] ) <= 0:
                print('Difference between last sample and first in line search is 0')
                break


        # searching along larger alpha values
        if lsiter == 0:
            self.alpha_max = self.alpha_max / self.beta     ## update max step length
        
        elif lsiter > 2:
            self.alpha_max = self.alpha_max * self.beta**(lsiter-2)     ## update max step length

        ## in case of exhausting numer of iterations stop (TODO: probably not necessary)
        elif lsiter >= self.maxlsiter:
            print('Error - line search ... EXITING')
            self.alpha_max = self.alpha_max * self.beta
      
        ## save best result
        f_min = np.min(f)
        ix = np.argmin(f)
        alpha_min = alpha[ix]
        print('\tnum iter: %d' %(lsiter))
        
        return alpha_min, f_min


    ##
    #
    #
    #
    #
    def interpolateInitialization( self, badVols=None, X=None ):
        
        # load data
        if not X is None:
            X = self.X
        T = len(X)

        # create bad vol list if necessary
        if badVols is None:
            outlierMask_lowres = np.zeros([T])
            for tt in range(T):
                for spk in self.badSpokes:
                    if spk in self.spkOfVol[tt]:
                        outlierMask_lowres[tt] = 1
            badVols = np.where(outlierMask_lowres)[0]
        print(badVols)
        # define T axis and bad vols
        time = [tt for tt in range(T) if tt not in badVols ]
        print(time)
        goodX = [X[tt] for tt in range(T) if tt not in badVols ]

        # learn interpolation
        goodX = np.array(goodX).reshape( len(time), np.prod(self.imgDim) )
        interp_fcn = interp1d( time, goodX.T  ,kind='linear', bounds_error=False, fill_value='extrapolate')

        # interpolate
        new_X = interp_fcn(badVols).T

        # fill the gaps
        for tt in range(len(badVols)):
            X[badVols[tt]] = new_X[tt,:].reshape(self.imgDim)

        return X 


    #%%
    #
    #   Register the sequence of volumes. 
    #   This function saves the results in a temporary folder, applies the selected registration procedure and loads the results.
    #   multiple registration methods will be implemented. if none is specified no registration is done
    #       'standard' -> registration to single reference volume specified in referenceImgIX
    #       'sequential -> register every volume to the next volume
    #   TODO: allow for multiple registration methods
    #
    def registerVolumes(self, outputDir=None, regMethod=None, masksPath=None, referenceImgIX=None, endBaselineIX=0, sliceOverSampling=0 ):

        if outputDir is None:
            outputDir = self.tempDir

        ## save volume sequence
        graspVol_uint16, image4DFile, self.imageFiles = increaseDimVolume( np.array(self.X).transpose(1,2,3,0), outputDir +  'graspRecon', \
                        referenceNii=None, flipDim=1, sliceOverSampling=sliceOverSampling, newImgSize=self.imgDim, saveIndividualVols=True, rescaleBool=False  )
        
        # create mask
        masks = {}
        masks['body'] = outputDir + 'bodyMask.nii'
        bodyMask = np.ones(self.X[0].shape)
        niiRefFile = nib.load(  self.imageFiles[0] )
        mask_nii = nib.Nifti1Image( np.array(bodyMask, dtype=np.uint8), \
                            niiRefFile.affine, niiRefFile.header)
        mask_nii.to_filename(masks['body'] )

       
        if regMethod == 'standard':
            ## select reference image if unspecified
            # select a volume towards the end (more contrast should help)
            if referenceImgIX is None:
                referenceImgIX = np.round(self.numVol * 3/4)

            referenceImgIX = [ referenceImgIX for tt in range(self.numVol) ]
            initialTransform = self.registrationTFM

        elif regMethod == 'sequential':
            referenceImgIX = [ tt+1 for tt in range(self.numVol) ]
            referenceImgIX[self.numVol-1] = self.numVol-1
            initialTransform = None
        
        ## reset registration to original position (clean registration). 
        # Note that when this is called initialTransform should be None
        # This option is necessary to do sequential registration, otherwise the composition of registrations is a mess.
        # Parallel(n_jobs=self.num_cores)(delayed(applyRegistration) \
        #         (self.imageFiles[0], self.imageFiles[tt], self.inverseregistrationTFM[tt], self.imageFiles[tt], tt,  mode='elastix_dense'  ) for tt in range( self.numVol ) )
        
        ## save previous transform files
        tempTFM = [ outputDir + os.path.basename(self.registrationTFM[tt][:-4]) + '_prev.txt'  for tt in range(self.numVol) ]
        for tt in range(self.numVol):
            ii = 1
            while os.path.exists(tempTFM[tt]):
                tempTFM[tt] = outputDir + os.path.basename(self.registrationTFM[tt][:-4]) + '_prev%d.txt' %(ii)
                ii += 1
            copyfile( self.registrationTFM[tt], tempTFM[tt] )

        ## apply registration method
        if regMethod == 'RLTI':
            newTFM, finalGeneratedSequence, finalRegisteredSequence = runRLTIregistration( image4DFile, None, 'noFIDsignal', masksFolder=masksPath, output_folder=outputDir, \
                                referenceImgIX=referenceImgIX, badVolumes=[], endBaselineIX=endBaselineIX, num_stokes_per_image=34, num_cores=self.num_cores, tempFolder=None, maxIter=1,\
                                referenceNii=self.imageFiles[0] )
        else:
            newTFM, finalRegisteredSequence = \
                register_images2reference( outputDir, self.imageFiles, output_folder=outputDir, \
                                referenceIX=referenceImgIX, \
                                num_cores= self.num_cores, tempFolder=None, initialTransform=initialTransform, iterNum=self.regIterNum )

        ## compose with previous iterations
        for tt in range(self.numVol):
            insertComposedRegistrationParameterFile( newTFM[tt], self.registrationTFM[tt], composedTransformFile=tempTFM[tt])

        self.regIterNum += 1    # update iteration number
                        
        ## compute inverse registration
        for tt in range(self.numVol):
            self.inverseregistrationTFM[tt] = self.registrationTFM[tt][:-4] + '_inverse.txt'

        Parallel(n_jobs=self.num_cores)(delayed(computeInverseRegTransform) \
                (self.imageFiles[tt], self.registrationTFM[tt], self.inverseregistrationTFM[tt], tt) for tt in range( self.numVol ) )

        ## re-load volume sequence and flip the slice coordinate
        regData = nib.load(finalRegisteredSequence).get_fdata()
        self.X = [ np.fft.fftshift( regData[:,:,:,tt], axes=2) for tt in range(self.numVol) ]
        

    ##
    #
    #
    #
    def saveAndApplyRegistration( self, X, registrationTFM, filename=None, sliceOverSampling=0.2 ):
        
        if filename is None:
            tempfile = 'graspRecon'

        # temporarily save image
        graspVol_uint16, image4DFile, imageFiles = increaseDimVolume( np.array(X).transpose(1,2,3,0), self.tempDir +  tempfile, \
                    referenceNii=None, flipDim=1, sliceOverSampling=sliceOverSampling, newImgSize=self.imgDim, saveIndividualVols=True, rescaleBool=False, save4DFile=True  )

        # apply registration
        Parallel(n_jobs=self.num_cores)(delayed(applyRegistration) \
                (imageFiles[0], imageFiles[tt], registrationTFM[tt], imageFiles[tt][:-4] + 'registered', tt,  mode='elastix_dense'  ) for tt in range( self.numVol ) )

        # load image and return
        regData = nib.load(image4DFile).get_fdata()

        return [ np.fft.fftshift( regData[:,:,:,tt], axes=2) for tt in range(self.numVol) ]

    
    ##
    #
    #   This function runs GRASP reconstruction in a pyramidal way. 
    #   It starts creating a low resolution reconstruction and interpolates to a higher dimension until reaching the desired dimension.
    #
    #   TODO: so far it is implemented for 1D recon only
    #
    def pyramidInitialization( self, tempRes, lambdas ):

        numLev = len(tempRes)

        originalPhase = np.copy( self.spkPhase )
        originalSpkBinWith = np.copy( self.spkBinWidth )
        originalNumVol = len(self.X)

        self.defineSpokesOfVol( np.linspace(0,1,self.kDataDim[2]).reshape(self.kDataDim[2],1), [float(tempRes[0])] )
        self.X = self.iNUFFT()
        numVol = len(self.X)
        

        # for every pyramidal level
        for ll in range(numLev):
            
            print('Level %d. Temporal Resolution %0.6f. Num vol %d' %(ll, tempRes[ll], int( 1/float(tempRes[ll]) )) )

            # set-up phases of current level
            self.defineSpokesOfVol( np.linspace(0,1,self.kDataDim[2]).reshape(self.kDataDim[2],1), [float(tempRes[ll])] )

            # # update normalization
            # self.xNorm = np.prod(self.imgDim)*self.numVol

            # run GRASP recon
            self.solve_CG('GRASP', lamReg=lambdas[ll], maxIterCG=10, regMethod=None, endBaselineIX=0, referenceImgIX=None, masksPath=None) 
            
            # get sizes
            if ll+1 < numLev:
                numVol_next = int( 1/float(tempRes[ll+1]) )
            else:
                numVol_next  = originalNumVol

            # interpolate results
            func = interp1d( np.linspace(0,1,self.numVol), np.array(self.X).reshape( ( self.numVol, int(np.prod(self.imgDim)) ) ).T, kind='linear',fill_value='extrapolate')
            newInterpX = func( np.linspace(0,1,numVol_next))
            self.X = [ newInterpX[:,tt].reshape(self.imgDim) for tt in range(numVol_next) ]


        # return to desired size
        self.defineSpokesOfVol( originalPhase, originalSpkBinWith )
        # self.xNorm = np.prod(self.imgDim)*self.numVol 

