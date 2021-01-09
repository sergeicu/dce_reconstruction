#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 12 09:13:29 2018

@author: jaume
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

# math imports
import numpy as np
import scipy as sp
from scipy import signal
from sklearn.cluster import k_means
# import matplotlib.pyplot as plt
from joblib import Parallel, delayed
import nibabel as nib
import nrrd

# locally created imports
# sys.path.append( os.path.expanduser("~") +  '/Documents/Research/DCE/registration/')
sys.path.append('/home/ch199899/Documents/Research/DCE_all/registration/')
from sequential_registration import listValidFiles, loadListofFiles, saveVolumeToListOfFiles

#%% generateLTIImpulseResonse
#   
#    This function generates the impulse response of a system that can be represented with the sum of single real poles.
#    
#   The definition of the poles used in this function is taken from the paper:
#       B. Yilmaz, K. Bekiroglu, C. Lagoa, and M. Sznaier, “A Randomized Algorithm for Parsimonious Model Identification,” IEEE Trans. Automat. Contr., vol. 63, no. 2, pp. 532–539, 2017.
#
#   
def generateLTIImpulseResonse( poles, selSet, time, pthreats=2 ):
    
    P = poles.size
    N = time.size
    
    # for every pole
    parallelResults = Parallel(n_jobs=pthreats)( delayed(createSingleLTIResponse)(poles[pp], selSet[pp], time) for pp in range(P) )
    impulseResponse = list(range(P))
    LTIsystems  = list(range(P))
    for pp in range(P):

        impulseResponse[pp] = parallelResults[pp][1]
        LTIsystems[pp] = parallelResults[pp][0]

    return impulseResponse, LTIsystems

def createSingleLTIResponse( poles, selSet, time):

    N = time.size

    # compute scaling factor
    sign = np.random.randint(0,2,1)*2 -1
    if (selSet == 1 ) | (selSet == 2):
        phia = (1 - np.abs(poles)**(2*N)) / (1 - np.abs(poles)**(2))
        phip = (1 - poles**(2*N)) / (1 - poles**(2))
        Gamma = ( np.real(phip) - phia \
                    - np.real( poles**2*np.conj(poles)**(2*N)*phip ) \
                    + np.abs(poles)**(2*N+2)*phia ) \
                / ( 1 - np.abs(poles)**2 )

    if selSet == 3:
        alphap = sign*( 1- poles**2 ) / ( 1 - poles**(2*N+2) )
    elif (selSet == 1):
        alphap = sign*np.sqrt( 2*(np.real(phip**2) + phia**2) + 2*np.sqrt(2*Gamma*( np.abs(phip)**2 - phia**2 )) )**(-1)
    elif (selSet == 2):
        alphap = sign*np.sqrt( 2*(phia**2 - np.real(phip**2)) + 2*np.sqrt(2*Gamma*( np.abs(phip)**2 - phia**2 )) )**(-1)
    elif selSet == 4:
        alphap = sign*1.0

    # generate impulse response
    ltiPoles = [poles]
    if selSet == 1:   # complex with real coef
        alphap = 2*alphap
        ltiPoles = [poles, np.conj(poles)]
        zeros = [np.real(poles)]
    elif selSet == 2:   # complex with imag coef
        alphap =  2*np.imag(poles)*alphap
        ltiPoles = [poles, np.conj(poles)]
        zeros = [0.0]
    else:
        zeros = [0.0]

#### MUST CHANGE THE ZEROS FOR CASES 1 AND 2!!!!
    # generate LTI systems
    LTIsystems = signal.dlti( zeros, ltiPoles, alphap, dt=1 )
    
    # generate impulse response
    impulseResponse = LTIsystems.impulse(t=time)[1][0].reshape(N)
    
    if any( np.isnan(impulseResponse) ):
        print( 'Gato! %0.6f, j%0.6f' %( np.real(poles), np.imag(poles) ))

    return LTIsystems, impulseResponse

#%% random selection of atoms
# 
# 
# 
def samplePolesFromUnitDisk( Np, rhoMax=1.0-1e-6, phiMax=180 ):

    poles = np.zeros([Np],dtype=np.complex)
    # selSet = list(range(Np))

    # select set from which to sample
    if phiMax > 0:
        selSet = np.random.randint(1,5,Np)
    else:
        selSet = np.ones((Np,))*3  # exclude poles with imaginary part

    # for each pole
    for pp in range(Np-1):
        
        # sample pole for each case
        if selSet[pp] == 3:
            poles[pp] = np.random.uniform(0.,rhoMax,10)[0]
        elif (selSet[pp] == 1) | (selSet[pp] == 2):
            length = np.sqrt( np.random.uniform( 0.,rhoMax,10)[0] )
            angle = np.pi * np.random.uniform(-phiMax/180.,phiMax/180.,10)[0]
            poles[pp] = length * np.exp( 1j* angle )
        elif selSet[pp] == 4:   # TODO: double check that this is correct. Shouldn't this represent constant values?
            length = np.sqrt( np.random.uniform( 0.,rhoMax,10)[0] )
            angle = np.pi * np.random.uniform(-phiMax/180.,phiMax/180.,10)[0]#* np.random.uniform(-5/180.,5/180.,10)[0]
            poles[pp] = length * np.exp( 1j* angle )
        
    poles[Np-1] = 1
    selSet[Np-1] = 4
 
    return poles, selSet
    

#%% run LTI system identification method
#    This function runs the system identification method described in:
#        B. Yilmaz, K. Bekiroglu, C. Lagoa, and M. Sznaier, “A Randomized Algorithm for Parsimonious Model Identification,” IEEE Trans. Automat. Contr., vol. 63, no. 2, pp. 532–539, 2017.
#    It consists on convex optimization of the system using pre-defined atoms and a Frank-Wolfe algorithm.
#    The original paper proposes to use a random set of atoms. Here we use a fixed set
#
#
#
def runLTIsysID( y, P, sampleTime, tau, badVolumes=[], maxIter = 100, Tu=None, rhoMax=1.0-1e-6, phiMax=180.0 ):
    
    N = sampleTime.size
    stepEpsilon = 1e-6    
    pthreats = 20

    # select good volumes
    goodVolumes = np.ones([N])
    goodVolumes[badVolumes] = 0
    goodVolumes = np.where(goodVolumes)[0]
    
    # define input
    if Tu is None:
        Tu = np.eye(N)
    
    # compute initial xk estimate
    poleAtoms, selectionSet = samplePolesFromUnitDisk( 1,  rhoMax, phiMax )
    impulseResponse, LTIatoms = generateLTIImpulseResonse( poleAtoms, selectionSet, sampleTime, pthreats )
    xk = impulseResponse[-1]*tau #*0.0
    ck= [tau]
    
    # plt.figure()
    # ITERATE  !
    k = 0
    while True:
        
        ## ----- random atom selection step ----- ##
        poleAtoms, selectionSet = samplePolesFromUnitDisk( P, rhoMax )
        impulseResponse, LTIatoms = generateLTIImpulseResonse( poleAtoms, selectionSet, sampleTime )
        
        ## ----- precompute gradient ----- ##
        gradF = 2*Tu.T.dot(Tu.dot(xk) - y )
        
        ## ----- Minimize over all selected atoms ----- ##
        # for all candidate atoms compute inner product
        rp = Parallel(n_jobs=pthreats)( delayed(np.dot)(impulseResponse[pp][goodVolumes],gradF[goodVolumes]) for pp in range(P) )
        
        # select minimum projection
        minK = np.argmin(np.array(rp))

        ## ----- compute step length ----- ##
        a = Tu.dot(tau* impulseResponse[minK] - xk)
        alpha = max(min( -np.dot((Tu.dot(xk)[goodVolumes] - y[goodVolumes]), a[goodVolumes] ) / np.dot(a[goodVolumes],a[goodVolumes])  ,1) ,0)
        
        ## ----- update xk ----- ##
        x_prev = xk
        xk = xk + alpha*(tau* impulseResponse[minK] - xk)
        
        # plot 
#            #for ii in range(P):
    #     plt.clf()
    #     plt.plot(sampleTime,y)
    #     plt.plot(sampleTime,xk)
    #  #        plt.plot(impulseResponse[ii])
    #     plt.plot(sampleTime,gradF)
    #     plt.plot(sampleTime,impulseResponse[minK],'--')
    #     plt.legend(['y','xk','grad','imp'])
    #     plt.show()
    #     plt.pause(0.1)

        ## ----- update ck list ----- ##
        # ck[minK] += alpha
        
        ## ----- evaluate fit and convergence ----- ##
        fitErr = np.linalg.norm( (Tu.dot(xk) - y)[goodVolumes]  ,ord=2)
        stepSize = np.linalg.norm( (xk - x_prev)[goodVolumes], ord=2)
        # print( 'Iter: %d. Obj. fun.: %0.6f. Alpha: %0.6f  Step: %0.6f. Atom (%d): %0.6f + j%0.6f' %( k,  fitErr, alpha, stepSize,  minK , np.real(poleAtoms[minK]), np.imag(poleAtoms[minK]) ))
        
        if (k > maxIter) | (stepSize < stepEpsilon):
            # print( 'Gato Dominguez!')
            break

        
        # update counter
        k +=1 
     
    return ck, xk, LTIatoms
        
#%% run LTI system identification method for multiple vectors
#    
#   Runs the runLTIsysID function in parallel for multiple vectors.
#       Optionally, it can allow for pre-clustering to avoid computation at cost of accuracy
#
#   INPUT:
#       DATA:
#       - inputImgs - <Nx,Ny,Nz,T>double - these are the input images where to fit the LTI model
#       - masks - dict - dictionary containing an entry for each mask:
#                           - aorta - <Nx,Ny,Nz>{1/0} - aorta mask
#                           - tissue - <Nx,Ny,Nz>{1/0} - mask for the tissue to be fitted
#                   ASSUME NON-OVERLAPPING MASKS!!!
#       LTI FIT PARAMS:
#       - numAtoms  - int - number of atoms to sample in each iteration (default is 100)
#       - tau       - int - regularization parameter in LTI FIT  default(100)
#       - badVolumes  - list - list of time instances considered outliers (will be ignored during fitting)  (default [])
#       - maxIter     - int - maximum number of iterations to run LTI (default 1000)
#       - numClusters - int - number of K-means clusters to use for fitting (defaut 100)
#       DCE SPECIFIC:
#       - aif   - <1,T>double - artery input function, if defined this one will be used. (default None)
#       - endBaselineIX - int - time instance (in entry position) when the contrast agent is injected. Used to determine aif (default 0)
#
#   OUTPUT:
#       - genDataTissue - <Nx,Ny,Nz,T> double - LTI fit to the provided input data. Masked with the mask['tissue'].
#       - aif - <1,T>double - artery input function, if defined this one will be used. (default None)
#
def runLTISysID4DCE( inputImgs, masks, numAtoms=100, tau=100, badVolumes=[], maxIter=1000, numClusters=100, nthreads=20, aif=None, endBaselineIX=0, rhoMax=1.0-1e-6, phiMax=180.0  ):    

    ## params
    T = inputImgs.shape[-1]

    dictKeys = list(masks.keys())    # explicitly stating that is a list is necessary for python 3
    ## run recon on aorta and compute aif
    if aif is None:
        if any( 'aifMask' == s for s in dictKeys  ):

            print( 'Computing aif from aifMask data')
            # mask aorta
            maskIX = np.where(masks['aifMask'].ravel())[0]
            dataMtrx = inputImgs.reshape( np.prod(inputImgs.shape[:-1]),T )[maskIX,:] 
            
            print(badVolumes)
            badVolumes_AIF = []
            for bv in badVolumes:
                print(bv)
                if bv >= endBaselineIX:
                    badVolumes_AIF.append( bv-endBaselineIX )

            # run LTI reconstruction with 10 clusters
            genDataAorta, parresults, labels = groupLTISysID( dataMtrx[:,endBaselineIX:], numAtoms,  tau, \
                        badVolumes=badVolumes_AIF, maxIter=maxIter, numClusters=10, nthreads=nthreads, inSignal=None, rhoMax=rhoMax, phiMax=phiMax )

            # fill the missing AIF with the average value
            genDataAorta = np.concatenate((  np.tile(np.mean(dataMtrx[:,:endBaselineIX],axis=1),[endBaselineIX,1]).T ,genDataAorta ), axis=1)


            # compute mean signal
            # aif = np.mean( genDataAorta,axis=0 )

            # pick the largest signal in the aorta
            clmax = np.argmax( np.max(genDataAorta, axis=-1) )
            aif =  genDataAorta[clmax,:]

            dictKeys.remove('aifMask')

        else:
            print( 'No aorta mask found. Using impulse at %d as aif.' %(endBaselineIX))
            aif = np.zeros([1,T])
            aif[0,endBaselineIX] = 1.0

    ## run recon on all other masks
    ltiImg = np.zeros((np.prod(inputImgs.shape[:-1]), T))
    for key in dictKeys:

        print( 'Running ' + key + ' group from ' + ','.join(dictKeys))
        # mask tissue
        maskIX = np.where(masks[key].ravel())[0]
        dataMtrx = inputImgs.reshape( np.prod(inputImgs.shape[:-1]) ,T )[maskIX,:] 
        
        # run LTI reconstruction with 10 clusters
        if maskIX.size > 0:
            fittedLTI, parresultsTissue, labels = groupLTISysID( dataMtrx, numAtoms,  tau, badVolumes=badVolumes, \
                            maxIter=maxIter, numClusters=numClusters, nthreads=nthreads, inSignal=aif)
            
            ltiImg[maskIX,:] = fittedLTI


    # reshape results
    ltiImg = ltiImg.reshape(inputImgs.shape) 

    return ltiImg, aif



#%% run LTI system identification method for multiple vectors
#    
#   Runs the runLTIsysID function in parallel for multiple vectors.
#       Optionally, it can allow for pre-clustering to avoid computation at cost of accuracy
#
#
def groupLTISysID( dataMtrx, numAtoms, tau, badVolumes=[], maxIter=100, numClusters=100, nthreads=20, inSignal=None, rhoMax=1.0-1e-6, phiMax=180.0  ):    

    # normalize all voxels
    T = dataMtrx.shape[-1]
    minY = np.min(dataMtrx,axis=-1)
    normY = ( np.max(dataMtrx,axis=-1) - minY )
    normY[normY==0] = 1
    dataMtrx = (dataMtrx-np.tile(minY,[T,1]).T)/np.tile(normY,[T,1]).T

    if not inSignal is None:
        inSignal = ( inSignal - np.min(minY) ) / np.max( normY )

    # apply clustering
    if dataMtrx.shape[0] > numClusters:
        
        sizeGroup = int(10*numClusters)
        numGroups = int(np.ceil(dataMtrx.shape[0] / float(sizeGroup)))

        print('Clustering Data into %d groups of size %d' %(numGroups, sizeGroup))
        if dataMtrx.shape[0] >= numGroups*numClusters:

            # divide data into sets for easier computation of clusters (for computational reasons)
            # sizeGroup = np.ceil( dataMtrx.shape[0] / float(numGroups) ) 
            groupmeanSignal = np.zeros((0,T))
            groupLabels = list(range(numGroups))
            groupIX = list(range(numGroups))
            for gg in range(numGroups):
                groupIX[gg] = np.arange( gg*sizeGroup, min( dataMtrx.shape[0], (gg+1)*sizeGroup ) , dtype=np.int )
                if groupIX[gg].size > numClusters:
                    tmpMeanSignal, groupLabels[gg] = clusterData( dataMtrx[groupIX[gg],:], numClusters=numClusters, nthreads=nthreads )   
                else:
                    tmpMeanSignal = dataMtrx[groupIX[gg],:]
                    groupLabels[gg] = np.arange(groupIX[gg].size)

                groupmeanSignal = np.concatenate( (groupmeanSignal,tmpMeanSignal), axis=0 )
        else:
            numGroups = 1
            groupIX = [np.arange(dataMtrx.shape[0])]
            groupmeanSignal = dataMtrx

        # now compute the cluster of clusters
        meanSignal, labelClusters = clusterData( groupmeanSignal, numClusters=numClusters, nthreads=nthreads )

        if dataMtrx.shape[0] > numGroups*numClusters:

            # reassign labels 
            groupFinalLabels = list(range(numGroups))
            for gg in range(numGroups):
                labelIX = np.arange( gg*numClusters, min( groupmeanSignal.shape[0], (gg+1)*numClusters ) , dtype=np.int )
                groupFinalLabels[gg] = np.zeros((groupIX[gg].size,))
                for cl in range(numClusters):
                    for cl2 in np.where(labelClusters[labelIX] == cl)[0]:
                        groupFinalLabels[gg][ groupLabels[gg] == cl2 ] = cl

            labels = np.zeros((dataMtrx.shape[0],))
            # reassign labels 
            for gg in range(numGroups):
                labels[ groupIX[gg] ] = groupFinalLabels[gg]
                
        else:
            labels = labelClusters

    else:
        meanSignal = dataMtrx
        labels = np.arange(dataMtrx.shape[0])

    Nclust = meanSignal.shape[0]


    # for every curve, fit
    print( 'Fitting LTI model' )
    parresults = Parallel(n_jobs=nthreads)( delayed(normLTIit) \
            ( meanSignal[cl,:], numAtoms,  tau, badVolumes, maxIter, inSignal=inSignal, rhoMax=rhoMax, phiMax=phiMax ) for cl in range(Nclust))
        
    # recontruct data matrix
    outData = np.zeros(dataMtrx.shape)
    for cl in range(Nclust):
        outData[ labels == cl ,:] = parresults[cl]

    # de-normalize voxels
    outData = np.tile(normY,[T,1]).T * outData + np.tile(minY,[T,1]).T

    return outData, parresults, labels


#%% normalized run LTI
#
#       This function takes an input signal, normalizes it and fits the LTI model
#
#
def normLTIit( y, numAtoms, tau, badvols=[], maxIter=100, inSignal=None, rhoMax=1.0-1e-6, phiMax=180.0  ):

    T = y.size

    minY = np.min(y)
    normY = ( np.max(y) - np.min(y) )
    if normY==0:
        normY = 1

    y = (y-minY)/normY
    
    if inSignal is None:
        Tu = np.eye(T)
    else:
        inSignal = (inSignal - minY)/normY 
        Tu = sp.linalg.toeplitz(inSignal, np.zeros((T,)))

    genData = runLTIsysID( y, numAtoms, np.arange(T), tau, badvols, maxIter=maxIter, Tu=Tu, rhoMax=rhoMax, phiMax=phiMax )[1]

    out = Tu.dot(genData)
    
    out = out*normY + minY
    
    return out


#%% 
#       This function computes a K-means clustering of the provided data to perform the fitting more efficiently
#
#
#
def clusterData( dataMtrx, numClusters=100, nthreads=20 ):

    centroid, labels = k_means(dataMtrx, n_clusters=numClusters, n_jobs=nthreads, precompute_distances=False)[:2]

    return centroid, labels



## MAIN FUNCTION
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='')
    parser.add_argument('-i', '--input_folder', required=True,
                        help='input folder to the nrrd sequence')
    parser.add_argument('-o', '--output_folder', required=True,
                        help='output folder')
    parser.add_argument('-m', '--mask', required=True,
                        help='path to mask')
    parser.add_argument('-tau', '--tau', default=100.0,
                        help='Regularization parameter Tau (default 100.0)')
    parser.add_argument('-numAtoms', '--numAtoms', default=100,
                        help='Number of atoms sampled in each iteration (default 100)')
    parser.add_argument('-tag', '--tag', required=True,
                        help='name tag of the nrrd sequence')
    parser.add_argument('-notag', '--notag', default='',
                        help='Non desired tag in the FSL files to load. (e.g. to exclude 4D volume)')
    parser.add_argument('-p', '--cores', default=20,
                        help='Number of cores to use (default 20)')
    args = parser.parse_args()

    #%% Set params
    input_folder = args.input_folder
    output_folder = args.output_folder
    mask_path = args.mask
    tau = float(args.tau)
    num_cores = int(args.cores)
    numAtoms = int(args.numAtoms)

    #%% load mask file (resample if needed)
    mask = nib.load( mask_path ).get_fdata()
    kidneyIX = np.where( mask.ravel() )[0]
    Nx_m,Ny_m,Nz_m = mask.shape


    #%% load ordered sequence of images
    imageFiles = listValidFiles( input_folder, goodTag=args.tag, badTag=args.notag, posIX=[7] )
    dataMtrx, imgHeaders = loadListofFiles( input_folder, imageFiles, mask )

    #%% run rLTI
    print( 'Running LTI. This may take a while...')
    maskIX = np.where(mask.ravel())[0]
    Nvx = maskIX.size
    T = dataMtrx.shape[-1]
    print( num_cores)
    parresults = Parallel(n_jobs=num_cores)( delayed(runLTIsysID) \
            ( dataMtrx[vx,:]/np.max(np.abs(dataMtrx[vx,:])), numAtoms, np.arange(T), tau ) for vx in range(Nvx))

    #%% parse results
    Nx, Ny, Nz = mask.shape
    recImg = np.zeros( (Nx*Ny*Nz, T) )
    for vx in range(Nvx):
        recImg[kidneyIX[vx],:] = parresults[vx][1]*np.max(dataMtrx[vx,:])
    recImg[np.isnan(recImg)] = -1

    #%% save results to images
    generatedImages = list(range(T))
    for tt in range(T):
        splitRoot = os.path.splitext(os.path.basename( imageFiles[tt] ))
        generatedImages[tt] = output_folder + splitRoot[0] + '_LTIgenerated' + splitRoot[-1]
        
    saveVolumeToListOfFiles( recImg, generatedImages, imgHeaders, mask )
    
    # create 4D volume
    joinImageSequence( generatedImages, output_folder + '/LTIgeneratedJoined.nrrd' )