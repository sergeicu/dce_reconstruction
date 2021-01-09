#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 11 11:00:56 2019

@author: ch199899
"""
import sys
import os

import numpy as np 
import scipy.io as spio
import matplotlib.pyplot as plt 
import nibabel as nib 
import hdf5storage



sys.path.insert(0, os.path.expanduser("~") + '/Documents/Research/DCE/DCEmodel/')
from DCE_tools import signal_concentration
sys.path.insert(0, os.path.expanduser("~") + '/Documents/Research/DCE/image_reconstruction/')
from DCE_tools import detectCorruptedSpokes



#%% define metrics
def takeImagePatch(volume, width, pos):
    xx,yy = np.meshgrid(pos[0]  + np.arange(-width/2,width/2),  pos[1] + np.arange(-width/2,width/2))
    return volume[xx,yy]

def compute_psnr( volume, sigMask, bkgMask ):
    xx,yy = np.where( bkgMask )
    mse = np.mean(  volume[ xx,yy]**2  )
    xx,yy = np.where( sigMask )
    return  np.max( volume[ xx,yy]**2 ) / mse 
    
def compute_snr( volume, sigMask, bkgMask ):
    xx,yy = np.where( bkgMask )
    noise = np.std(  volume[ xx,yy]  )
    xx,yy = np.where( sigMask )
    S = np.mean(  volume[ xx,yy]  )
    return 0.655 * S / noise 

def compute_coefVar( volume, sigMask, bkgMask ):
    xx,yy = np.where( sigMask )
    return np.std(  volume[ xx,yy ]  ) / np.mean(  volume[ xx,yy ]  )
    
def compute_tsnr( volume, sigMask ):
    xx,yy = np.where( sigMask )
    patchVol = volume[xx,yy,:]
    meanSig = np.mean( patchVol )
    stdSig = np.std( patchVol )
    return  meanSig / stdSig


 
#%% 
def plotAreaWithMetrics( volumes, width, pos, bkgPos, titleSTR, startContrast=0, mask=None, widthBkg=None, time=None, clim=None):

    zz = np.array(pos[2],dtype=np.int)
    
    Nx,Ny,Nz,T = volumes[0].shape
    
    if mask is None:
        xx, yy = np.meshgrid(pos[0]  + np.arange(-width/2,width/2),  pos[1] + np.arange(-width/2,width/2))
    else:
        xx, yy = np.where(mask[:,:,zz])

    if widthBkg is None:
        widthBkg = width
    
    if time is None:
        time = range(T)

    if clim is None:
        clim = (0,np.max(volumes[0].ravel()))
        

    
#    bkgX, bkgY = np.meshgrid( np.arange(0,50), np.arange(0,50) )
    bkgX, bkgY = np.meshgrid( bkgPos[0]  + np.arange(-widthBkg/2,widthBkg/2),  bkgPos[1] + np.arange(-widthBkg/2,widthBkg/2))
    bkgMask = np.zeros((Nx,Ny,Nz))
    bkgMask[np.array(bkgX,dtype=np.int),np.array(bkgY,dtype=np.int) ,np.array(zz,dtype=np.int)] = 1
    
    psnr = np.zeros((T,len(volumes)))
    SNR = np.zeros((T,len(volumes)))
    tsnr = np.zeros((T,len(volumes)))
    coefVar = np.zeros((T,len(volumes)))
    plt.figure()
    for tt in time:
        plt.clf()
        for ss in range(len(volumes)):
            plt.subplot(1,len(volumes),ss+1)
            voxelMask = np.zeros((Nx,Ny,Nz))
            voxelMask[np.array(xx,dtype=np.int),np.array(yy,dtype=np.int),np.array(zz,dtype=np.int)] = 1
            plt.imshow( volumes[ss][:,:,zz,tt].T, cmap='gray' ,clim=clim)
            plt.contour(voxelMask[:,:,zz].T,colors='r')
            plt.contour(bkgMask[:,:,zz].T,colors='g')
            plt.xticks([])
            plt.yticks([])
            
            psnr[tt,ss] = compute_psnr( volumes[ss][:,:,zz,tt]*(bkgMask+voxelMask)[:,:,zz], voxelMask[:,:,zz], bkgMask[:,:,zz] )
            SNR[tt,ss] =  compute_snr( volumes[ss][:,:,zz,tt]*((bkgMask+voxelMask)[:,:,zz]), voxelMask[:,:,zz], bkgMask[:,:,zz] )
            tsnr[tt,ss] =  compute_tsnr( volumes[ss][:,:,zz,:], voxelMask[:,:,zz] )
            coefVar[tt,ss] = compute_coefVar( volumes[ss][:,:,zz,tt]*((bkgMask+voxelMask)[:,:,zz]), voxelMask[:,:,zz], bkgMask[:,:,zz] )
            metricResultsSTR = titleSTR[ss] + ' ' + str(tt) + '\n PSNR: %0.2f\n SNR: %0.2f\n tSNR: %0.2f\n coefVar: %0.2f' %(psnr[tt,ss], SNR[tt,ss], tsnr[tt,ss], coefVar[tt,ss])

            plt.title(metricResultsSTR)
        plt.pause(0.2)
        
    print('Signal-To-Noise results:')
    print('\nSNR:')
    plt.figure()
    plt.plot(SNR)
    plt.title('SNR')
    plt.legend( titleSTR  )
    print('Difference before contrast %0.3f +/-%0.3f' %(np.mean( np.diff( SNR[:startContrast,:] ,axis=1) ), np.std( np.diff( SNR[startContrast:,:] ,axis=1) )))
    print('Difference after contrast %0.3f +/-%0.3f' %(np.mean( np.diff( SNR[startContrast:,:] ,axis=1) ), np.std( np.diff( SNR[startContrast:,:] ,axis=1) )))
    print('Difference overall %0.3f +/-%0.3f' %(np.mean( np.diff( SNR ,axis=1) ), np.std( np.diff( SNR[startContrast:,:] ,axis=1) )))
    
    print('\nPSNR:')
    plt.figure()
    plt.plot(psnr)
    plt.title('PSNR')
    plt.legend( titleSTR  )
    print( 'Difference before contrast %0.3f +/-%0.3f' %(np.mean( np.diff( psnr[:startContrast,:] ,axis=1) ), np.std( np.diff( psnr[startContrast:,:] ,axis=1) )))
    print( 'Difference after contrast %0.3f +/-%0.3f' %(np.mean( np.diff( psnr[startContrast:,:] ,axis=1) ), np.std( np.diff( SNR[startContrast:,:] ,axis=1) )))
    print( 'Difference overall %0.3f +/-%0.3f' %(np.mean( np.diff( psnr ,axis=1) ), np.std( np.diff( psnr[startContrast:,:] ,axis=1) )))
    
    print( '\nCV:')
    plt.figure()
    plt.plot(coefVar)
    plt.title('CV')
    plt.legend( titleSTR  )
    print( 'Difference before contrast %0.3f +/-%0.3f' %(np.mean( np.diff( coefVar[:startContrast,:] ,axis=1) ), np.std( np.diff( coefVar[startContrast:,:] ,axis=1) )))
    print( 'Difference after contrast %0.3f +/-%0.3f' %(np.mean( np.diff( coefVar[startContrast:,:] ,axis=1) ), np.std( np.diff( coefVar[startContrast:,:] ,axis=1) )))
    print( 'Difference overall %0.3f +/-%0.3f' %(np.mean( np.diff( coefVar ,axis=1) ), np.std( np.diff( coefVar[startContrast:,:] ,axis=1) )))
    
    return psnr, SNR, tsnr, coefVar





#%%------- LOAD DATA ---------#
plt.close('all')

#%% set subject info
subjects = [
{'patientName':'NAME_SURNAME', 'mrn':1111, 'endBaselineIX':7}
]


spk_per_vol=34
shadedAreacolor = 'green'
shadedStartcolor = 'xkcd:turquoise'


MFT = []


for subj in subjects:
    
    #%% load
    baseFolder = '/fileserver/projects6/jcollfont/DCE/Processed/'
    graspData = baseFolder + 'Subj' + str(subj['mrn']) + '/' + subj['patientName'] + '/common_processed/GRASP_reconstruction_autoBadCoils/subj'+str(subj['mrn'])+'_recon_GRASP_N448_S34_lambda0.012500_vol4D.nii'
    elstxData = baseFolder + 'Subj' + str(subj['mrn']) + '/' + subj['patientName'] + '/common_processed/GRASP_registered_elastix/result.0.nii.gz'
    rltiData = baseFolder + 'Subj' + str(subj['mrn']) + '/' + subj['patientName'] + '/common_processed/GRASP_registered_RLTI_kmeansAll/registered_body_iter3.nii.gz'
    refVolData = baseFolder + 'Subj' + str(subj['mrn']) + '/' + subj['patientName'] + '/common_processed/GRASP_registered_RLTI_kmeansAll/registered_init_body.nii.gz'
    ltiVolData = baseFolder + 'Subj' + str(subj['mrn']) + '/' + subj['patientName'] + '/common_processed/GRASP_registered_RLTI_kmeansAll/generated_body_iter2.nii.gz'
    ltiInterpVolData = baseFolder + 'Subj' + str(subj['mrn']) + '/' + subj['patientName'] + '/common_processed/GRASP_registered_RLTI_kmeansAll/registered_interp_body_iter3.nii.gz'
    maskPaths = os.path.dirname(graspData) + '/'
    FIDpath = baseFolder + 'Subj' + str(subj['mrn']) + '/' + subj['patientName'] + '/common_processed/' + 'FID_folder_fast/subj' + str(subj['mrn']) + '_FID_all_avrg_signal.mat'
    
    
    dataPaths = {'GRASP':graspData,'REFVOL':refVolData,'ELSTX':elstxData, 'RLTI':rltiData}
    
#    regMethods = ['GRASP','REFVOL','ELSTX','RLTI']
    regMethods = ['GRASP','RLTI']
    
    #%% load image sequences
    data = {}
    for rm in regMethods:
        data[rm] = nib.load(dataPaths[rm]).get_fdata()
        
    imgSize = data['GRASP'].shape
    
    #%% laod masks
    masks = {}
    for mm in ['aorta','rightKidney','leftKidney']:
        maskName = maskPaths + mm + 'Mask.nii'
        if os.path.exists( maskName ):
            masks[mm] = nib.load(maskName).get_fdata()
    
    #%% detect outliers
    outlierMaskVolIX, corrFID_thr = detectCorruptedSpokes(FIDpath, spk_per_vol, subj['endBaselineIX']*1.5)
    
    MFT.append( 1- len(outlierMaskVolIX) / imgSize[-1] )
    
   