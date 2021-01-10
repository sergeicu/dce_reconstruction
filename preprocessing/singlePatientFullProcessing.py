#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 15:45:18 2017

@author: ch194093
"""

import numpy as np
import os 
import sys
from subprocess import call
#sys.path.insert(0, '/home/ch194093/Desktop/kidneydcemri/kidneySegmentation/')
from initProcessFuncs import dicomHeaderInfoExtractor,initProcess
import pandas as pd
import nibabel as nib
import re
#sys.path.insert(0, '/home/ch194093/Desktop/kidneydcemri/preprocessing-visualization')
import funcs

import pydicom
import shutil


from plotTools import plotMIPVideo
    
#sys.path.insert(0, '/fileserver/external/rdynaspro4/abd/deepLearningSegmentation')
#import singlePatientSegMethods
os.environ["CUDA_DEVICE_ORDER"] = "PCI_BUS_ID"   # see issue #152
os.environ["CUDA_VISIBLE_DEVICES"] = "0"  ##P4000=='0', P6000=='1'

############################# Flags ####################
downloadDataEnabled=1;
preprocessingEnabled=1;
visualizationEnabled=0;
createMIPs=1;
kidneyMaskEnabled=0;
aortaMaskEnabled=0;
########################################################

############################# Kidney segmentation params ####################
kidneyMaskParams={};
kidneyMaskParams['TestSetNum']=1;
kidneyMaskParams['tpUsed']=50;
kidneyMaskParams['tDim']=kidneyMaskParams['tpUsed'];
kidneyMaskParams['PcUsed']=0;
kidneyMaskParams['networkToUse']='Unet';
kidneyMaskParams['visualizeResults']=visualizationEnabled;  
kidneyMaskParams['visSlider']=visualizationEnabled;
########################################################







# Get the Inputs: patient "name", "MRN" and "date"
def inputInfo():
    personName = input('Enter patient Name (First_Last):')
    personMRN = input('Enter patient MRN:')
    personDate = input('Enter Acquisition Date (yyyymmdd):')
    print(personName+'  '+personMRN+'  '+personDate+'  ');

    confirmVar = input('Is the above information correct? (y/n):')
    if confirmVar=='n':
        personName, personMRN, personDate=inputInfo();
    return personName, personMRN, personDate


#### Read subject list and add the entered information to the list

fileAddress='/fileserver/external/rdynaspro4/abd/MRUcommon/subjectList_MRU.xls';
subjectInfo=pd.read_excel(fileAddress, sheet_name=0)

if downloadDataEnabled:    
    personName, personMRN, personDate=inputInfo();        
    
    
    if personName not in subjectInfo:
        subjectInfo[personName]=None;
      
        
    subjectInfo.loc[0,personName]=personDate;
    subjectInfo.loc[1,personName]=personMRN;
    
    
    ###### Download data from server 

    patientName=personName;
    date=str(subjectInfo._get_value(0,patientName))
    MRN=str(subjectInfo._get_value(1,patientName))
    outPutDir='/fileserver/commondataraw/MRUoriginal/'+patientName;
    os.system('mkdir '+outPutDir);
    #oscommand='/fileserver/abd/MRUcommon/retrieve.sh '+MRN+' '+date+' '+outPutDir+' MR';

    oscommand='/fileserver/external/rdynaspro4/abd/MRUcommon/retrieve2.sh '+MRN+' '+date+' '+outPutDir+' MR';

    os.system(oscommand);
    firstSubFolder=os.listdir(outPutDir);
    print(firstSubFolder)

    if len(firstSubFolder[0])>8:
        firstSubFolder1=firstSubFolder[0][0:7]+'?'+firstSubFolder[0][8:]
        os.system('mv '+outPutDir+'/'+firstSubFolder1+' '+outPutDir+'/'+personMRN);

else:
    personName, personMRN, personDate=inputInfo();  
    
    patientName=personName;
        
    if personName not in subjectInfo:
        subjectInfo[personName]=None;
      
        
    subjectInfo.loc[0,personName]=personDate;
    subjectInfo.loc[1,personName]=personMRN;
    
    
print('Output folder is:')    
print('/fileserver/commondataraw/MRUoriginal/'+patientName)
input('Execution paused. Delete bad folders')
    
###### find data path and sequence prefixes

fileAddress='/fileserver/external/rdynaspro4/abd/MRUcommon/subjectDicomInfo.xls';
subjectInfo2=pd.read_excel(fileAddress, sheet_name=0)

if preprocessingEnabled:
    
    diir0='/fileserver/commondataraw/MRUoriginal/'+patientName;
    
    diir1=sorted(os.listdir(diir0+'/'+os.listdir(diir0)[0]))

    diir2=os.listdir(diir0+'/'+os.listdir(diir0)[0]+'/'+diir1[-1]) #change back sila
    #diir2=os.listdir(diir0+'/'+os.listdir(diir0)[0]+'/'+diir1[0])
    diir2=[diir2[i] for i in range(len(diir2)) if diir2[i].split("_")[0].isdigit()];
#    diir2writeTime=[os.path.getmtime(diir0+'/'+os.listdir(diir0)[0]+'/'+diir1[-1]+'/'+diir2[i]) for i in range(len(diir2))];
#    argSortTime=np.argsort(diir2writeTime);
#    diir22=[diir2[argSortTime[i]] for i in range(len(diir2))];
#    diir2=np.copy(diir22);
        
    argSortStartTime=np.argsort([int(f.split("_")[0]) for f in diir2]);
    diir22=[diir2[argSortStartTime[i]] for i in range(len(diir2))];
    diir2=np.copy(diir22);
    
    diir0+'/'+os.listdir(diir0)[0]+'/'+diir1[-1]
    
    diir3=[f for f in diir2 if "_T" in f and ((len(f) - 1 - f[::-1].index("T_"))>len(f)-4) and \
               str.isdigit(f[(len(f) - 1 - f[::-1].index("T_"))+1])]
    diir3_2=[f for f in diir2 if "_TT" in f and f.index("_TT")>len(f)-12]
    flag_TT=0;
    if len(diir3)<=len(diir3_2):
        diir3=diir3_2.copy();
        flag_TT=1;
    
    
    if (len(diir3) == 0) & (len(diir3_2) == 0):      ## TODO: @JCF introduced this if to account for the recons from VITA. The folder structure changes
        if any([ f.find('DYNAMIC_GRASP_ONUR') > -1 for f in diir22 ] ):    # the sequence contains the right name
            ## TODO: @JCF is there a folder with the name matching the GRASP sequence in VITA?
            
            graspFolders = np.where([ diir22[f]  == '%d_OBL_COR_STARVIBE_6min_DYNAMIC_GRASP_ONUR' %(int(diir22[f].split('_')[0])) for f in range(len(diir22)) ])[0]
            if graspFolders.size == 0:  # in case the data is from VIDA
                graspFolders = np.where([ diir22[f]  == '%d_OBL_COR_STARVIBE_DYNAMIC_GRASP_ONUR' %(int(diir22[f].split('_')[0])) for f in range(len(diir22)) ])[0]
            dicomFodlerPath = diir0+'/'+os.listdir(diir0)[0]+'/'+diir1[-1]+ '/'+ diir22[graspFolders[0]] + '/'
            dicomFiles = os.listdir( dicomFodlerPath )   ## assuming there is one and only one folder matching
            
            
            ## retrieve the acquisition time of all slices
            sliceLocations = []
            acqTime = []
            for dd in dicomFiles:
                dicomData = pydicom.read_file( dicomFodlerPath + dd ,force=True)
                
                sliceLocations.append(dicomData.SliceLocation)
                time =dicomData.AcquisitionTime
                acqTime.append(int(time.split(".")[0][0:2])*3600+int(time.split(".")[0][2:4])*60+int(time.split(".")[0][4:6])+(float(time)-int(float(time))))
                
            # create new folder structure
            diir3 = []
            unique_ACQtimes = np.unique( acqTime )
            timeRes = np.diff(unique_ACQtimes)
            for uACQ in range(len(unique_ACQtimes)):
                
                patPath = diir0+'/'+os.listdir(diir0)[0]+'/'+diir1[-1] + '/'
                
                newDir = '%d_COR_OBL_STARVIBE_DYNAMIC_dyn_3sec_TT=%s' %( len(diir2) +1 + uACQ,  unique_ACQtimes[uACQ] )
                if uACQ > 0:
                    if timeRes[uACQ-1] > 10:
                        newDir =  '%d_COR_OBL_STARVIBE_DYNAMIC_dyn_12sec_TT=%s' %( len(diir2) +1 + uACQ,  unique_ACQtimes[uACQ] )
                        
                diir3.append( newDir )
                if not os.path.exists( newDir ):
                    os.makedirs( patPath + diir3[-1] )
                    
                for dd in range(len(dicomFiles)):
                    if acqTime[dd] == unique_ACQtimes[uACQ]:
                        shutil.copyfile( dicomFodlerPath + dicomFiles[dd], patPath + diir3[-1] +'/'+ dicomFiles[dd] )
                 
            flag_TT=0;
                
        
    diir4=[len(f) for f in diir3];
    uniqElements= list(set(diir4))
    countRepetition=[diir4.count(i) for i in uniqElements];
    indicesOfElementsWith1rep = [diir4.index(uniqElements[i]) for i, x in enumerate(countRepetition) if x == 1 \
                                 and ((uniqElements[i]+1 not in uniqElements) and (uniqElements[i]-1 not in uniqElements))];
    wrongElements=[diir3[i] for i in indicesOfElementsWith1rep]
    if len(uniqElements)>1:
        for i2 in range(len(indicesOfElementsWith1rep)):
            diir3.remove(wrongElements[i2]);
                
        
#    diir3=sorted(diir3);
    underLineNdx1=[m.start() for m in re.finditer('_', diir3[0])]
    underLineNdx2=[m.start() for m in re.finditer('_', diir3[-1])]
    seqPer1=diir3[0][underLineNdx1[0]:underLineNdx1[-1]+2]
    seqPer2=diir3[-1][underLineNdx2[0]:underLineNdx2[-1]+2]
    if flag_TT==1:
        seqPer1=diir3[0][underLineNdx1[0]:underLineNdx1[-1]+3]
        seqPer2=diir3[-1][underLineNdx2[0]:underLineNdx2[-1]+3]
    if seqPer1==seqPer2: #added by sila for singel sequence cases
        seqs=seqPer1
    else:
        seqs=[seqPer1,seqPer2]
            
    
            
    
    subjectInfo2.loc[patientName] = 0;
    subjectInfo2.loc[[patientName],['patPath']]=diir0+'/'+os.listdir(diir0)[0]+'/'+diir1[-1]+'/' #sila change
    #subjectInfo2.loc[[patientName],['patPath']]=diir0+'/'+os.listdir(diir0)[0]+'/'+diir1[0]+'/'
    if len(list(set(seqs)))==2:
        subjectInfo2.loc[[patientName],['seqPrefix1']]=seqs[0]
        subjectInfo2.loc[[patientName],['seqPrefix2']]=seqs[1]
        subjectInfo2.loc[[patientName],['numSeq']]=len(list(set(seqs)))
    else:
        subjectInfo2.loc[[patientName],['seqPrefix1']]=seqs
        subjectInfo2.loc[[patientName],['seqPrefix2']]=None;
        subjectInfo2.loc[[patientName],['numSeq']]=1
    
    ###################### run init process
    correctTime="False";fixFOV="True";fixIntJump="False";scaleFactor=1;
    
    patPath=subjectInfo2['patPath'][patientName]
    seqPrefix1=subjectInfo2['seqPrefix1'][patientName]
    seqPrefix2=subjectInfo2['seqPrefix2'][patientName]
    numSeq=subjectInfo2['numSeq'][patientName]
    initProcess(patientName,patPath,seqPrefix1,seqPrefix2,correctTime,fixFOV,fixIntJump,scaleFactor,numSeq)
    
    ############################ save 4D nii
    scanRecDir='/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/';
    scanRecDirSubs=os.listdir(scanRecDir);
    # scanRecDirSubs=scanRecDirSubs[5:];
    scanRecDirSubs=scanRecDirSubs;
    matchingDirs = [s for s in scanRecDirSubs if patientName in s];
    matchingDirs=sorted(matchingDirs)
    images={}
    for s in matchingDirs:
        subjectDirNii=[i for i in os.listdir(scanRecDir+s) if os.path.isfile(os.path.join(scanRecDir+s,i)) and 'reconSCAN_T' in i];
        print(s)
        
    #     print(str(len(subjectDirNii)))
        for k in range(len(subjectDirNii)):   
             s1='reconSCAN_T'+str(k)+'.nii'
             print(s1)
             img0 = nib.load(scanRecDir+s+'/'+s1)
             data = img0.get_data()
             data1=np.expand_dims(data,axis=3)
             if s1[11]==str(0):
                img=data1;
             else:
                img=np.concatenate((img, data1), axis=3);
#        for s1 in subjectDirNii:
#    #         print(s1[-5])
#            img0 = nib.load(scanRecDir+s+'/'+s1)
#            data = img0.get_data()
#            data1=np.expand_dims(data,axis=3)
#            if s1[11]==str(0):
#                img=data1;
#                print('yes')
#            else:
#                img=np.concatenate((img, data1), axis=3);
        images[s]=img;    
        print(img.shape)    
        img4D = nib.Nifti1Image(img, img0.affine)
        nib.save(img4D,scanRecDir+s+'/reconSCAN_4D.nii')
    
    ######################### Extract Dicom info
    
    
    x=dicomHeaderInfoExtractor(patientName,patPath, seqPrefix1, seqPrefix2,numSeq);
    subjectInfo2['timeRes'][patientName]=x['timeRes'];
    subjectInfo2['StationName'][patientName]=x['StationName'];
    subjectInfo2['ModelName'][patientName]=x['ManufacturerModelName'];
    subjectInfo2['PatientID'][patientName]=x['PatientID'];
    subjectInfo2['PatientAge'][patientName]=x['PatientAge'];
    subjectInfo2['seq1EndTime'][patientName]=x['seq1EndTime'];
    subjectInfo2['seq2StartTime'][patientName]=x['seq2StartTime'];
    subjectInfo2['seq1VolNum'][patientName]=x['seq1VolNum'];
    subjectInfo2['seq2VolNum'][patientName]=x['seq2VolNum'];
    
    # TODO: @JCF. commented out the next two lines. They breaak the code and might be leftovers from previous hacks to deal with pandas and XLS file
    subjectInfo[patientName][6]=x['seq1VolNum'];
    subjectInfo[patientName][7]=x['seq2VolNum'];
    
    ####################################### save excel files
    writer = pd.ExcelWriter('/fileserver/external/rdynaspro4/abd/MRUcommon/subjectDicomInfo.xls')
    subjectInfo2.to_excel(writer,'Sheet1')
    writer.save();
    
    
    writer = pd.ExcelWriter('/fileserver/external/rdynaspro4/abd/MRUcommon/subjectList_MRU.xls')
    subjectInfo.to_excel(writer,'Sheet1')
    writer.save();
    if numSeq==2:
        oscommand='chmod -R 777 '+scanRecDir+patientName+'_seq1';
        os.system(oscommand);
        oscommand='chmod -R 777 '+scanRecDir+patientName+'_seq2';
        os.system(oscommand);
    else:
        oscommand='chmod -R 777 '+scanRecDir+patientName;
        os.system(oscommand);
############################## Visualization sequence1
if visualizationEnabled:
    
    
#    patientNameNum=124
    import matplotlib.pyplot as plt
    from matplotlib.widgets import Slider
    if preprocessingEnabled:
        im=images[matchingDirs[0]];  
        imgPath = '/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/'+patientName+'_seq1/reconSCAN_4D.nii'
    else:
        #if subjectInfo2['numSeq'][patientName]==2:
        if subjectInfo2['numSeq'][patientName]==2:  
            imgPath = '/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/'+patientName+'_seq1/reconSCAN_4D.nii'
        else:
            imgPath = '/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/'+patientName+'/reconSCAN_4D.nii'
        
        img= nib.load(imgPath)
        im=img.get_data()  

    maxRange=im.max()-20
    fig1 = plt.figure(1, figsize=(6,6))
    ax1 = fig1.add_subplot(111)
    fig1.subplots_adjust(left=0.25, bottom=0.25)
    fig1.suptitle('Sequence1')
    s0 = 2;
    t0 = 0;
    
    # im = max0 * np.random.random((10,10))
    ax1.imshow(im[:,:,s0,t0].T,origin='upper',aspect='auto',cmap='gray',clim=(0, maxRange))
    # fig.colorbar(im1)
    
    axcolor = 'lightgoldenrodyellow'
    axmin = fig1.add_axes([0.25, 0.1, 0.65, 0.03])
    axmax  = fig1.add_axes([0.25, 0.15, 0.65, 0.03])
    
    st1 = Slider(axmin, 'Time (0:'+str(im.shape[3])+')', 0, im.shape[3], valinit=t0,valfmt='%d')
    ss1 = Slider(axmax, 'Slice', 0, im.shape[2], valinit=s0,valfmt='%d')
    
    def update(val):
        timeNdx=st1.val;
        sliceNdx=ss1.val;       
        ax1.imshow(im[:,:,int(sliceNdx),int(timeNdx)].T,origin='upper',aspect='auto',cmap='gray',clim=(0, maxRange))
    #     im1.set_clim([smin.val,smax.val])
        fig.canvas.draw()
        
    st1.on_changed(update)
    ss1.on_changed(update)
    
    plt.show()
    
    
    ################################Visualization seq2
    if subjectInfo2['numSeq'][patientName]==2:
        if preprocessingEnabled:
            im2=images[matchingDirs[1]];
            imgPath = '/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/'+patientName+'_seq2/reconSCAN_4D.nii'
        else:
            img2= nib.load('/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/'+patientName+'_seq2/reconSCAN_4D.nii')
            im2=img2.get_data()  
                        
        
        maxRange=im2.max()-20
        fig = plt.figure(2, figsize=(6,6))
        ax = fig.add_subplot(111)
        fig.subplots_adjust(left=0.25, bottom=0.25)
        fig.suptitle('Sequence2')
        s0 = 2;
        t0 = 0;
        
        # im = max0 * np.random.random((10,10))
        ax.imshow(im2[:,:,s0,t0].T,origin='upper',aspect='auto',cmap='gray',clim=(0, maxRange))
        # fig.colorbar(im1)
        
        axcolor = 'lightgoldenrodyellow'
        axmin = fig.add_axes([0.25, 0.1, 0.65, 0.03])
        axmax  = fig.add_axes([0.25, 0.15, 0.65, 0.03])
        
        st = Slider(axmin, 'Time (0:'+str(im2.shape[3])+')', 0, im2.shape[3], valinit=t0,valfmt='%d')
        ss = Slider(axmax, 'Slice', 0, im2.shape[2], valinit=s0,valfmt='%d')
        
        def update(val):
            timeNdx=st.val;
            sliceNdx=ss.val;       
            ax.imshow(im2[:,:,int(sliceNdx),int(timeNdx)].T,origin='upper',aspect='auto',cmap='gray',clim=(0, maxRange))
        #     im1.set_clim([smin.val,smax.val])
            fig.canvas.draw()
            
        st.on_changed(update)
        ss.on_changed(update)
        
        plt.show();


if createMIPs:
    
    
    if subjectInfo2['numSeq'][patientName]==2:  
        imgPath = ['/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/'+patientName+'_seq1/reconSCAN_4D.nii','/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/'+patientName+'_seq2/reconSCAN_4D.nii']
    else:
        imgPath = ['/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/'+patientName+'/reconSCAN_4D.nii']
    
    for imP in imgPath:
        
        img= nib.load(imP)
        im=img.get_data()  
        
        plotMIPVideo( [im], frameRateX=4, titles=None, maxVal=None, vidName= os.path.dirname(imP) +  '/MIP_video.mp4' )
    
    


if 0:
### volume intensity drift check
    imVec=np.mean(np.reshape(im,(im.shape[0]*im.shape[1]*im.shape[2],im.shape[3])),0);
    im2Vec=np.mean(np.reshape(im2,(im2.shape[0]*im2.shape[1]*im2.shape[2],im2.shape[3])),0);
    fid=np.concatenate((imVec,im2Vec),axis=0);
    plt.figure();plt.plot(fid,'-*');
    #im2=+(imVec[-1]-im2Vec[0]);
    ############# find the baseline for aorta
    #from collections import Counter
    baseline=funcs.baselineFinder(im);
#    print(baseline)
#    aortaPotentialTimesIM=im[75:150,:,0:20,0:50]; # keep first 15 dataPoints
#    #aortaPotentialTimesIM=im[:,:,:,0:1];
#    x=(aortaPotentialTimesIM>.8*np.max(aortaPotentialTimesIM)).nonzero()    
#    
#    x=(np.max(aortaPotentialTimesIM,axis=3)-np.min(aortaPotentialTimesIM,axis=3)>.6*np.max(aortaPotentialTimesIM)).nonzero()    
#    #b = Counter(x[2]);
#    #mostOccInZ=b.most_common(1)[0][0];
#    
#    medianOfValsInXaxisNdx=(abs(x[0]-np.median(x[0]))<10).nonzero()[0];
#    medianOfValsInYaxisNdx=(abs(x[1]-np.median(x[1]))<10).nonzero()[0];
#    medianOfValsInZaxisNdx=(abs(x[2]-np.median(x[2]))<10).nonzero()[0];
#    commonXYconstraint=np.intersect1d(medianOfValsInXaxisNdx,medianOfValsInYaxisNdx);
#    commonXYconstraint=np.intersect1d(medianOfValsInZaxisNdx,commonXYconstraint)
#    allAortaPotentials=aortaPotentialTimesIM[x[0][commonXYconstraint],x[1][commonXYconstraint],x[2][commonXYconstraint],:]
#    
#    v,timeNdx = allAortaPotentials.max(1),allAortaPotentials.argmax(1)
#    maxA=np.median(timeNdx)
#    y=np.mean(allAortaPotentials[:,0:int(maxA)+1],0);
#    plt.figure();plt.plot(y3,'-*');
#    plt.figure();plt.plot(np.mean(allAortaPotentials,0),'-*');
##    fix min and max to [50 to 400] range
#    y2=(y*400)/y.max();y3=y2-y2.min()+50;
#    baseLine=(y3<0.5*(max(y3)-min(y3))).nonzero()[0][-1];
#    
#    print(baseLine)
#################  

#plt.close('all')



if kidneyMaskEnabled:
    box,image,zDimOrig,subjectInfo=singlePatientSegMethods.singlePatientDeepCassDet(personName,kidneyMaskParams);
    im,mask=singlePatientSegMethods.singlePatientDeepCassSeg(personName,box, image, kidneyMaskParams,zDimOrig,subjectInfo); 
    
    
if aortaMaskEnabled:
    call(['python2','AortaSeg_commandLine.py','-pn',personName, '-ns', str(subjectInfo2['numSeq'][personName])])
