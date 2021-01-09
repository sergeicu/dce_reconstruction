#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 08:05:48 2019

@author: ch199899
"""

import sys
import os
import numpy as np 
import scipy.stats as spst
import scipy.io as spio
import nibabel as nib
import nrrd
import pandas as pd
import ast
from skimage.morphology import closing, opening
from copy import deepcopy
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as manimation

sys.path.insert(0, os.path.expanduser("~") + '/Documents/Research/DCE/DCEmodel/')
from DCE_tools import signal_concentration

sys.path.insert(0, os.path.expanduser("~") + '/Documents/Research/DCE_all/model_fitting/')
from sourbronModelTools import estimateAIF




def plotImagesWithMaskVideo( data, regKeys, masks, timeIX, sliceIX, timeRes=[3.3], maxValIn=None, vidName=None ):
    
    if not vidName is None:
        FFMpegWriter = manimation.writers['ffmpeg']
        metadata = dict(title='Recon plots', artist='',
                        comment='')
        writer = FFMpegWriter(fps=15, metadata=metadata)
    
    T = -1
    timeIX = list(range(len(regKeys)))
    for rm in range(len(regKeys)):
        T = np.maximum(data[regKeys[rm]].shape[-1],T)
        timeIX[rm] = np.linspace(0,1,data[regKeys[rm]].shape[-1])
        
    timeAll = np.linspace(0,1,T)
    
    fig, ax = plt.subplots(1,len(regKeys))
    if len(regKeys) == 1:
        ax = [ax]
    
    with writer.saving(fig, vidName, 256):
        for tt in timeAll:
            
            plt.cla()
            for stale in range(1): # reduce frame rate
            
                for rm in range(len(regKeys)):
                    
                    imgTimeIX = np.argmin( np.abs( timeIX[rm] - tt  ) )
                    image = data[regKeys[rm]][:,:,sliceIX,imgTimeIX]
                    
                    if maxValIn is None:
                        maxVal = np.max(data[regKeys[rm]])
                    else:
                        maxVal = maxValIn[rm]
                    
                    
                    # plot image
                    ax[rm].imshow( image.T  \
                               ,clim=(0,maxVal), cmap='gray',alpha=1)
                    ax[rm].set_title( regKeys[rm])# + ' ' + str(tt))
                    
                    # plot mask contours
                    for km in masks.keys():
                        ax[rm].contour(masks[km][:,:,sliceIX].T, levels=[1],colors='r')
                    ax[rm].set_xticks( [] )
                    ax[rm].set_xticklabels( [] )
                    ax[rm].set_yticklabels( [] )
                    ax[rm].set_yticks( [] )
                
                if not vidName is None:
                    writer.grab_frame()
                
                
                
                
                
                
#%% 
subjects = [
{'patientName':None, 'mrn':1111111,'endBaselineIX':1},
]                
                
spk_per_vol=34
shadedAreacolor = 'green'
shadedStartcolor = 'xkcd:turquoise'            
                
#%%     
for subj in subjects:

    plt.close('all')

    try:
        #%% load
        baseFolder = '/fileserver/projects6/jcollfont/DCE/Processed/'
        subj['patientName'] = os.listdir( baseFolder + 'Subj' + str(subj['mrn']) + '/' )[0]
    #    graspData = baseFolder + 'Subj' + str(subj['mrn']) + '/' + subj['patientName'] + '/common_processed/RLTIrecon/graspRecon_4D.nii'
        graspData = baseFolder + 'Subj' + str(subj['mrn']) + '/' + subj['patientName'] + '/common_processed/GRASP_reconstruction_autoBadCoils/subj'+str(subj['mrn'])+'_recon_GRASP_N448_S34_lambda0.012500_vol4D.nii'
        elstxData = baseFolder + 'Subj' + str(subj['mrn']) + '/' + subj['patientName'] + '/common_processed/GRASP_registered_elastix/result.0.nii.gz'
        rltiData = baseFolder + 'Subj' + str(subj['mrn']) + '/' + subj['patientName'] + '/common_processed/GRASP_registered_RLTI_kmeansAll/registered_body_iter3.nii.gz'
        refVolData = baseFolder + 'Subj' + str(subj['mrn']) + '/' + subj['patientName'] + '/common_processed/GRASP_registered_RLTI_kmeansAll/registered_init_body.nii.gz'
        ltiVolData = baseFolder + 'Subj' + str(subj['mrn']) + '/' + subj['patientName'] + '/common_processed/GRASP_registered_RLTI_kmeansAll/generated_body_iter0.nii.gz'
        ltiInterpVolData = baseFolder + 'Subj' + str(subj['mrn']) + '/' + subj['patientName'] + '/common_processed/GRASP_registered_RLTI_kmeansAll/registered_interp_body_iter3.nii.gz'
        maskPaths = os.path.dirname(graspData) + '/'
        FIDpath = baseFolder + 'Subj' + str(subj['mrn']) + '/' + subj['patientName'] + '/common_processed/' + 'FID_folder_fast/subj' + str(subj['mrn']) + '_FID_all_avrg_signal.mat'
    
            #%% load XLS data
        xlsFile = '/fileserver/external/rdynaspro4/abd/MRUcommon/subjectDicomInfo.xlsx'
        xlsData = pd.read_excel(xlsFile)
        subjInfo = xlsData.loc[xlsData[xlsData.keys()[6]] == subj['patientName']]
        timeRes = np.array(subjInfo['timeRes'])[0]
        if not isinstance(timeRes,float):
            timeRes = ast.literal_eval(timeRes)
            scanData = '/fileserver/abd/GraspRecons/reconResultsSCAN/'+ subj['patientName']+ '_seq1/reconSCAN_4D.nii'
        else:
            scanData = '/fileserver/abd/GraspRecons/reconResultsSCAN/'+ subj['patientName']+ '/reconSCAN_4D.nii'
        
        dataPaths = {'GRASP':graspData,'REFVOL':refVolData,'ELSTX':elstxData, 'RLTI':rltiData, 'SCAN':scanData}
    
        regMethods = ['GRASP', 'SCAN']
    
    
        #%% load image sequences
        data = {}
        for rm in regMethods:
            
            if rm == 'SCAN':
                if not isinstance(timeRes,float):
                    data[rm] = nib.load('/fileserver/abd/GraspRecons/reconResultsSCAN/'+ subj['patientName']+ '_seq1/reconSCAN_4D.nii').get_fdata()
                    tempDat = nib.load('/fileserver/abd/GraspRecons/reconResultsSCAN/'+ subj['patientName']+ '_seq2/reconSCAN_4D.nii').get_fdata()
                    data[rm] = np.concatenate((data[rm],tempDat), axis=3)
                    
            else:
                data[rm] = nib.load(dataPaths[rm]).get_fdata()
            
            
            if data[rm].shape[0] == 448:
                data[rm] = data[rm][ np.arange(448/4,448/4*3,dtype=np.int),:,:,: ]
                data[rm] = data[rm][:,np.arange(448/4,448/4*3,dtype=np.int),:,:]
    
        imgSize = data['GRASP'].shape
        
        masks = {}
        for mm in ['aorta','rightKidney','leftKidney']:
            maskName = maskPaths + mm + 'Mask.nii'
            if os.path.exists( maskName ):
                masks[mm] = nib.load(maskName).get_fdata()
                
                
        
    #    #%%
    #    plotImagesWithMaskVideo( data,['GRASP'], masks, range(imgSize[-1]), int(imgSize[-2]/2), timeRes=[3.3], maxVal=np.max(data['GRASP']) ,vidName='graspVideo_'+ str(subj['mrn']) +'.mp4')
    #    
    #    #%%
    #    imgSize = data['SCAN'].shape
    #    plotImagesWithMaskVideo( data,['SCAN'], masks, range(imgSize[-1]), int(imgSize[-2]/2), timeRes=timeRes, maxVal=np.max(data['SCAN']) ,vidName='scanVideo_'+ str(subj['mrn']) +'.mp4')
    #    
        #%%
        plotImagesWithMaskVideo( data,['SCAN','GRASP'], {}, [np.max(data['SCAN']), np.max(data['GRASP'])], int(imgSize[-2]/2), timeRes=[3.3], maxValIn=None ,vidName='scannerVSgrasp_'+ str(subj['mrn']) +'.mp4')
        
    except:
        print('Subject %d failed' %(subj['mrn']))    