#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul  6 16:53:46 2017

@author: ch194093
"""
import nibabel as nib
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import signal
from scipy.ndimage import zoom
#import SimpleITK as sitk
#from matplotlib.widgets import Slider, Button, RadioButtons


def readData(patientName,subjectInfo,reconMethod,genBoundBox):
    
    seqNum=subjectInfo['numSeq'][patientName];

    # reconstruction Method
    # scanner ---> 'SCAN'
    # grasp   ---> 'GRASP'
    

    if reconMethod=='GRASP':
        dataAddress0='/fileserver/external/rdynaspro4/abd/GraspRecons/reconResults'+reconMethod+'/method2/';
    elif reconMethod=='SCAN':
        dataAddress0='/fileserver/external/rdynaspro4/abd/GraspRecons/reconResults'+reconMethod+'/';
#        dataAddress0='/fileserver/abd/GraspRecons/reconResults'+reconMethod+'/';


    if seqNum==1 or reconMethod=='GRASP':
        dataAddress=dataAddress0+patientName+'/recon'+reconMethod+'_4D.nii';
        img= nib.load(dataAddress)
        im=img.get_data()  
        
        maskAddress='/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/'+patientName+'/';
#        maskAddress='/fileserver/abd/GraspRecons/reconResultsSCAN/'+patientName+'/';
        if os.path.isfile(maskAddress+'aortaMask.nii'):
            am1= nib.load(maskAddress+'aortaMask.nii');am=am1.get_data()    
            if os.path.isfile(maskAddress+'leftKidneyMask.nii.gz'):
                lkm1= nib.load(maskAddress+'leftKidneyMask.nii.gz');lkm=2*lkm1.get_data()  
            else:
                lkm=np.zeros(np.shape(am));
                
            if os.path.isfile(maskAddress+'rightKidneyMask.nii.gz'):
                rkm1= nib.load(maskAddress+'rightKidneyMask.nii.gz');rkm=rkm1.get_data()  
            else:
                rkm=np.zeros(np.shape(am));
                
#            if np.sum(rkm[:])==0:
#                rkm=lkm;
#            if np.sum(lkm[:])==0:
#                lkm=rkm;
        else:
            lkm=0;rkm=0;am=0;
        
        
    elif seqNum==2:
        dataAddress=dataAddress0+patientName+'_seq1/recon'+reconMethod+'_4D.nii';  
        img1= nib.load(dataAddress)
        im1=img1.get_data()  
        dataAddress=dataAddress0+patientName+'_seq2/recon'+reconMethod+'_4D.nii';  
        img= nib.load(dataAddress)
        im2=img.get_data()  
    #     im=np.concatenate((im1, im2), axis=3);

        x=subjectInfo['timeRes'][patientName];
        seq1tres=float(x.split("[")[1].split(",")[0]);
        seq2tres=float(x.split(",")[1].split("]")[0]);
        if seq2tres>seq1tres:
            #resample second to first
            num2=seq2tres/seq1tres;
            im3=zoom(im2,(1,1,1,num2),order=0);
#            num=int(np.round(im2.shape[3]*seq2tres/seq1tres))
#            im3=signal.resample(im2, num, t=None, axis=3)
            im=np.concatenate((im1, im3), axis=3);
        else:
            im=np.concatenate((im1, im2), axis=3);
            
        maskAddress='/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/'+patientName+'_seq1/';
#        maskAddress='/fileserver/abd/GraspRecons/reconResultsSCAN/'+patientName+'_seq1/';
        if os.path.isfile(maskAddress+'aortaMask.nii.gz'):
            am1= nib.load(maskAddress+'aortaMask.nii.gz');am=am1.get_data()    
            if os.path.isfile(maskAddress+'leftKidneyMask.nii.gz'):
                lkm1= nib.load(maskAddress+'leftKidneyMask.nii.gz');lkm=2*lkm1.get_data()  
            else:
                lkm=np.zeros(np.shape(am));
                
            if os.path.isfile(maskAddress+'rightKidneyMask.nii.gz'):
                rkm1= nib.load(maskAddress+'rightKidneyMask.nii.gz');rkm=rkm1.get_data()  
            else:
                rkm=np.zeros(np.shape(am));
                
#            if np.sum(rkm[:])==0:
#                rkm=np.zeros(np.shape(lkm));
#            if np.sum(lkm[:])==0:
#                lkm=np.zeros(np.shape(rkm));   
        else:
            lkm=0;rkm=0;am=0;
    # dataAddress='/common/abd/GraspRecons/reconResults'+reconMethod+'/'+patientName+seqNum+'/recon'+reconMethod+'_4D.nii';
    # dataAddress='/common/abd/GraspRecons/reconResults'+reconMethod+'/'+patientName+'/recon'+reconMethod+'_4D.nii';
    # im=np.concatenate((im1, im2), axis=3); 
    boxes=[];
    if genBoundBox:
        aL=np.nonzero(lkm==2);
        aR=np.nonzero(rkm==1);
#        a1=np.zeros((np.shape(KM)))
#        a1[min(a[0]):max(a[0]),min(a[1]):max(a[1]),min(a[2]):max(a[2])]=4;
        # bounding box center + width and height
        if aL[0].size!=0:
            boxL=np.array([int((min(aL[0])+max(aL[0]))/2),int((min(aL[1])+max(aL[1]))/2),int((min(aL[2])+max(aL[2]))/2),\
              (max(aL[0])-min(aL[0])),(max(aL[1])-min(aL[1])),(max(aL[2])-min(aL[2]))])
        else:
            boxL=np.zeros((6,));
            
        if aR[0].size!=0:
            boxR=np.array([int((min(aR[0])+max(aR[0]))/2),int((min(aR[1])+max(aR[1]))/2),int((min(aR[2])+max(aR[2]))/2),\
              (max(aR[0])-min(aR[0])),(max(aR[1])-min(aR[1])),(max(aR[2])-min(aR[2]))])
        else:
            boxR=np.zeros((6,));
        
        boxes=np.vstack([np.array(boxR),np.array(boxL)]);
        
#    KM=np.logical_or(rkm, lkm);
    KM=rkm+lkm;
    im=(im/np.amax(im))*100;
    
    return im, KM, boxes


def writeMasks(patientName,subjectInfo,reconMethod,Masks2Save,overwrite):
    
    seqNum=subjectInfo['numSeq'][patientName];

    # reconstruction Method
    # scanner ---> 'SCAN'
    # grasp   ---> 'GRASP'
    

    if reconMethod=='GRASP':
        dataAddress0='/fileserver/external/rdynaspro4/abd/GraspRecons/reconResults'+reconMethod+'/method2/';
    elif reconMethod=='SCAN':
        dataAddress0='/fileserver/external/rdynaspro4/abd/GraspRecons/reconResults'+reconMethod+'/';
#        dataAddress0='/fileserver/abd/GraspRecons/reconResults'+reconMethod+'/';


    if seqNum==1 or reconMethod=='GRASP':
        dataAddress=dataAddress0+patientName+'/recon'+reconMethod+'_T0.nii';
        img= nib.load(dataAddress)
#        im=img.get_data()  
        
        maskAddress='/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/'+patientName+'/';
#        maskAddress='/fileserver/abd/GraspRecons/reconResultsSCAN/'+patientName+'/';
        if os.path.isfile(maskAddress+'leftKidneyMask.nii') and not overwrite:
#            nib.save(Masks2Save,maskAddress+'leftKidneyMask.nii');
#            nib.save(Masks2Save,maskAddress+'rightKidneyMask.nii');
#            am1= nib.save(maskAddress+'aortaMask.nii');am=am1.get_data()      
            print('Mask is already existant!');
        else:
            Masks2SaveR=Masks2Save['R'];Masks2SaveL=Masks2Save['L'];
            Masks2Save1R = nib.Nifti1Image(Masks2SaveR, img.affine)
            Masks2Save1L = nib.Nifti1Image(Masks2SaveL, img.affine)
            nib.save(Masks2Save1L,maskAddress+'leftKidneyMask_automatic.nii.gz');
            nib.save(Masks2Save1R,maskAddress+'rightKidneyMask_automatic.nii.gz');
#            am1= nib.save(maskAddress+'aortaMask.nii');am=am1.get_data()    
        
        
    elif seqNum==2:
        dataAddress=dataAddress0+patientName+'_seq1/recon'+reconMethod+'_T0.nii';  
        img1= nib.load(dataAddress)
#        im1=img1.get_data()  
#        dataAddress=dataAddress0+patientName+'_seq2/recon'+reconMethod+'_4D.nii';  
#        img= nib.load(dataAddress)
#        im2=img.get_data()  
    #     im=np.concatenate((im1, im2), axis=3);

#        x=subjectInfo['timeRes'][patientName];
#        seq1tres=float(x.split("[")[1].split(",")[0]);
#        seq2tres=float(x.split(",")[1].split("]")[0]);
#        if seq2tres>seq1tres:
#            #resample second to first
#            num=int(np.round(im2.shape[3]*seq2tres/seq1tres))
#            im3=signal.resample(im2, num, t=None, axis=3)
#            im=np.concatenate((im1, im3), axis=3);
#        else:
#            im=np.concatenate((im1, im2), axis=3);
            
        maskAddress='/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/'+patientName+'_seq1/';
#        maskAddress='/fileserver/abd/GraspRecons/reconResultsSCAN/'+patientName+'_seq1/';
        if os.path.isfile(maskAddress+'leftKidneyMask.nii') and not overwrite:
#            nib.save(Masks2Save,maskAddress+'leftKidneyMask.nii');
#            nib.save(Masks2Save,maskAddress+'rightKidneyMask.nii');
#            am1= nib.save(maskAddress+'aortaMask.nii');am=am1.get_data()      
            print('Mask is already existant!');
        else:
            Masks2SaveR=Masks2Save['R'];Masks2SaveL=Masks2Save['L'];
            Masks2Save1R = nib.Nifti1Image(Masks2SaveR, img1.affine)
            Masks2Save1L = nib.Nifti1Image(Masks2SaveL, img1.affine)
            nib.save(Masks2Save1L,maskAddress+'leftKidneyMask_automatic.nii.gz');
            nib.save(Masks2Save1R,maskAddress+'rightKidneyMask_automatic.nii.gz');
            oscommand='chmod -R 777 '+maskAddress;
            os.system(oscommand);
    # dataAddress='/common/abd/GraspRecons/reconResults'+reconMethod+'/'+patientName+seqNum+'/recon'+reconMethod+'_4D.nii';
    # dataAddress='/common/abd/GraspRecons/reconResults'+reconMethod+'/'+patientName+'/recon'+reconMethod+'_4D.nii';
    # im=np.concatenate((im1, im2), axis=3); 
#    KM=np.logical_or(rkm, lkm);
#    im=(im/np.amax(im))*100;
    
    return