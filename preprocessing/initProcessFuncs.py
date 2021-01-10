import numpy as np
import os 
import glob
import subprocess
from decimal import Decimal,getcontext
from collections import Counter
#import checkInput as chec
import sys

#def initProcess(patPath, seqPrefix1, seqPrefix2="",correctTime="False",fixFOV="False",fixIntJump="False",scaleFactor=1,numSeq=2):

def initProcess(patientName,patPath, seqPrefix1, seqPrefix2,correctTime,fixFOV,fixIntJump,scaleFactor,numSeq):
    #     if scaleFactor<0:
    # sys.exit("Please check the value of the scaling factor.")
    
    
############################################# if num of sequence is 1    
    if numSeq==1:
        if not patPath.endswith("/"):
            patPath=patPath+"/"
        files=glob.glob(patPath+"*"+seqPrefix1+"*")
        files=sorted(files, key=lambda name: int(name[len(patPath):].split('_')[0]));
        if len(files)==0:
            sys.exit("No directories with the given sequence prefix have been found.")
        os.system("mkdir "+patPath+"DCE_MRI_corrected/")
        for f in files:
            os.system("cp -r "+f+" "+patPath+"DCE_MRI_corrected/"+f.split("/")[-1])


        rescaleRealToShort="/fileserver/external/rdynaspro4/abd/MRUcommon/rescaleRealToShort"
        

        #Creating Directories for the new cropped and registered directories
        if os.path.isdir(patPath+"DCE_MRI"):
            os.system("rm -r " +patPath+"DCE_MRI")
        os.system("mkdir " +patPath+"DCE_MRI")
        if os.path.isdir(patPath+"DCE_MRI_NRRD"):
            os.system("rm -r " +patPath+"DCE_MRI_NRRD")
        os.system("mkdir " +patPath+"DCE_MRI_NRRD")


        #Copy DICOM to DCE_MRI directoy
        print("\nCopying DICOM to DCE_MRI...")
        sequence1=[];sequence1_0=[];
        for filename in files:
            sequence1_0.append(filename);
#            for c in filename:
#                if c=="(":
#                    filename=filename.replace(c,'\(')
#                if c==")"   :
#                    filename=filename.replace(c,'\)')
            filename=filename.replace('(','\(')
            filename=filename.replace(')','\)')
            #print(filename)
            sequence1.append(filename)
            
            os.system("cp -r "+filename+" "+patPath+"DCE_MRI")
        print("Successfully copied first sequence to DCE_MRI.")


        #Converting DICOM to NRRD
        print("\nConverting DICOM to NRRD...")
        for f in sequence1:
#            print(f)
            os.system("rm "+patPath+"DCE_MRI/"+f.split("/")[-1]+"/.~1~")
            os.system("mkdir "+patPath+"DCE_MRI_NRRD/"+f.split("/")[-1])
            os.system("crlConvertDICOMMOSAIC -d "+patPath+"DCE_MRI/"+f.split("/")[-1]+" -p "+patPath+"DCE_MRI_NRRD/"+f.split("/")[-1]+"/")
#            print("crlConvertDICOMMOSAIC -d "+patPath+"DCE_MRI/"+f.split("/")[-1]+" -p "+patPath+"DCE_MRI_NRRD/"+f.split("/")[-1]+"/")            

        print("Successfully converted first sequence from DICOM to NRRD.")

        
        #Converting DICOM to NIFTY
        print("\nConverting NRRD to NIFTY...")
        #tNdx=range(len(sequence1))
        tNdex=-1;
        for f in sequence1:
    #         print(f)
    #         print(f.split("/")[-1])
            tNdex=tNdex+1;
            nrddFileName=os.listdir(patPath+"DCE_MRI_NRRD/"+sequence1_0[tNdex].split("/")[-1]);
            nrddFileName1=nrddFileName[0];
#            for c in nrddFileName1:
#                if c=="(":
#                    nrddFileName1=nrddFileName1.replace(c,'\(')
#                if c==")":
#                    nrddFileName1=nrddFileName1.replace(c,'\)')
            nrddFileName1=nrddFileName1.replace('(','\(')  
            nrddFileName1=nrddFileName1.replace(')','\)')
            nifPath="/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/"+patientName+"/";
            os.system("mkdir "+nifPath)
            os.system("crlConvertBetweenFileFormats -in "+patPath+"DCE_MRI_NRRD/"+f.split("/")[-1]+"/"+nrddFileName1+" -out "+nifPath+"reconSCAN_T"+str(tNdex)+".nii")
            #print("crlConvertBetweenFileFormats -in "+patPath+"DCE_MRI_NRRD/"+f.split("/")[-1]+"/"+nrddFileName1+" -out "+nifPath+"reconSCAN_T"+str(tNdex)+".nii")

        print("Successfully converted first sequence from NRRD to Nii.")      
        
        
############################################# if num of sequence is 2        
    elif numSeq==2:
        #checkInput=chec.checkInput
        #patPath=checkInput(patPath, seqPrefix1, seqPrefix2)
        if not patPath.endswith("/"):
            patPath=patPath+"/"
        files=glob.glob(patPath+"*"+seqPrefix1+"*")
        files=sorted(files, key=lambda name: int(name[len(patPath):].split('_')[0]));
        files2=glob.glob(patPath+"*"+seqPrefix2+"*")
        files2=sorted(files2, key=lambda name: int(name[len(patPath):].split('_')[0]));
        rescaleRealToShort="/fileserver/external/rdynaspro4/abd/MRUcommon/rescaleRealToShort"

        ###########Creating Directories for the new cropped and registered directories

        if os.path.isdir(patPath+"DCE_MRI"):
            os.system("rm -r " +patPath+"DCE_MRI")
        os.system("mkdir " +patPath+"DCE_MRI")
        if os.path.isdir(patPath+"DCE_MRI_corrected"):
            os.system("rm -r " +patPath+"DCE_MRI_corrected")
        os.system("mkdir " +patPath+"DCE_MRI_corrected")
        if os.path.isdir(patPath+"DCE_MRI_NRRD"):
            os.system("rm -r " +patPath+"DCE_MRI_NRRD")
        os.system("mkdir " +patPath+"DCE_MRI_NRRD")
        if os.path.isdir(patPath+"DCE_MRI_NRRD_resampled"):
            os.system("rm -r " +patPath+"DCE_MRI_NRRD_resampled")
        os.system("mkdir " +patPath+"DCE_MRI_NRRD_resampled")
        if os.path.isdir(patPath+"DCE_MRI_NRRD_scaled"):
            os.system("rm -r " +patPath+"DCE_MRI_NRRD_scaled")
        os.system("mkdir " +patPath+"DCE_MRI_NRRD_scaled")
        print("Necessary directories have been created.")

        ###########Copy DICOM to DCE_MRI directoy

        print("\nCopying DICOM to DCE_MRI...")
        sequence1=[];sequence1_0=[];
        for filename in files:
            sequence1_0.append(filename);
            #sequence1.append(filename)
#            for c in filename:
#                if c=="(":
#                    filename=filename.replace(c,'\(')
#                if c==")"   :
#                    filename=filename.replace(c,'\)')
            filename=filename.replace("(","\(")
            filename=filename.replace(")","\)")
            #print(filename)
            sequence1.append(filename)
            #print("cp -r "+filename+" "+patPath+"DCE_MRI")
            os.system("cp -r "+filename+" "+patPath+"DCE_MRI")
        print("Successfully copied first sequence to DCE_MRI.")

        
        
        sequence2=[];sequence2_0=[];
        for filename in files2:
            sequence2_0.append(filename);
            #sequence2.append(filename)
#            if "(" in filename:
#               filename=filename.replace('(','\(')
#            if ")" in filename:
#               filename=filename.replace(')','\)')
#            print(filename)
            filename=filename.replace("(","\(")
            filename=filename.replace(")","\)")
            sequence2.append(filename)
            os.system("cp -r "+filename+" "+patPath+"DCE_MRI")
        print("Successfully copied both sequences to DCE_MRI.")

        ############  Converting DICOM to NRRD

        print("\nConverting DICOM to NRRD...")
        for f in sequence1:
          #  print(f)
            os.system("mkdir "+patPath+"DCE_MRI_NRRD/"+f.split("/")[-1])
            os.system("crlConvertDICOMMOSAIC -d "+patPath+"DCE_MRI/"+f.split("/")[-1]+" -p "+patPath+"DCE_MRI_NRRD/"+f.split("/")[-1]+"/")
            #print("crlConvertDICOMMOSAIC -d "+patPath+"DCE_MRI/"+f.split("/")[-1]+" -p "+patPath+"DCE_MRI_NRRD/"+f.split("/")[-1]+"/")            

        print("Successfully converted first sequence from DICOM to NRRD.")
        for f in sequence2:
#            print(f)
            os.system("mkdir "+patPath+"DCE_MRI_NRRD/"+f.split("/")[-1])
            os.system("crlConvertDICOMMOSAIC -d "+patPath+"DCE_MRI/"+f.split("/")[-1]+" -p "+patPath+"DCE_MRI_NRRD/"+f.split("/")[-1]+"/")
        print("Successfully converted both sequences from DICOM to NRRD.")

        
        ###########Converting NRRD to NIFTY
        print("\nConverting NRRD to NIFTY...")
        tNdex=-1;
        for f in sequence1:
            tNdex=tNdex+1;
            nrddFileName=os.listdir(patPath+"DCE_MRI_NRRD/"+sequence1_0[tNdex].split("/")[-1]);
            nrddFileName1=nrddFileName[0];
#            for c in nrddFileName1:
#                if c=="(":
#                    nrddFileName1=nrddFileName1.replace(c,'\(')
#                if c==")":
#                    nrddFileName1=nrddFileName1.replace(c,'\)')
            nrddFileName1=nrddFileName1.replace('(','\(')
            nrddFileName1=nrddFileName1.replace(')','\)')
            nifPath="/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/"+patientName+"_seq1/";
            os.system("mkdir "+nifPath)
            os.system("crlConvertBetweenFileFormats -in "+patPath+"DCE_MRI_NRRD/"+f.split("/")[-1]+"/"+nrddFileName1+" -out "+nifPath+"reconSCAN_T"+str(tNdex)+".nii")



        print("Successfully converted sequence 1 from NRRD to Nii.")        
        
        
        
        #Copy NRRD files from the first sequence to DCE_MRI_corrected        

        #if scaleFactor==1:
        #
         #   print("\nCopying NRRD files from the first sequence to DCE_MRI_corrected...")
          #  for filename in glob.glob(patPath+seqPrefix1+"*"):
           #     os.system("cp -r "+filename+" "+patPath+"DCE_MRI_corrected")

            #for file in glob.glob(patPath+seqPrefix1+"*"):
             #   k=0
              #  for f in os.listdir(patPath+"DCE_MRI_corrected/"+file.split("/")[-1]):
               #     os.system("mv "+patPath+"DCE_MRI_corrected/"+file.split("/")[-1]+"/"+f+" "+patPath+"DCE_MRI_corrected/"+file.split("/")[-1]+"/"+"crl_im%04d.dcm" %k)    
                #    k+=1
            #print("Successfully copied NRRD files from the first sequence to DCE_MRI_corrected.")

        #Scaling NRRDs

        print("\nSampling the first sequence...")
        for file in sequence1:
            dirName=file.split("/")[-1]            
            os.system("mkdir "+patPath+"DCE_MRI_NRRD_scaled/"+dirName)
            dirName2=dirName
#            for c in dirName2:
#                if c=="\\" :
#                    dirName2=dirName2.replace(c,'')
            dirName2=dirName2.replace('\\','')       
            baseName=os.listdir(patPath+"DCE_MRI_NRRD/"+dirName2)[0]
            baseName2=baseName
#            for c in baseName2:
#                if c=="(":
#                    baseName2=baseName2.replace(c,'\(')
#                if c==")":
#                    baseName2=baseName2.replace(c,'\)') 
            baseName2=baseName2.replace('(','\(')  
            baseName2=baseName2.replace(')','\)') 
            os.system("crlImageAddMultiplyAdd "+patPath+"DCE_MRI_NRRD/"+dirName+"/"+baseName2+" 0 "+str(scaleFactor)+" 0 "+patPath+"DCE_MRI_NRRD_scaled/"+dirName+"/"+baseName2)    
        print("Successfully sampled the first sequence");
        

        #Resampling of the second sequence

        print("\nResampling the second sequence...")
        dir0=sequence1[0];
#        for c in dir0:
#                if c=="\\" :
#                    dir0=dir0.replace(c,'')
        dir0=dir0.replace('\\','')
        os.chdir(dir0)
        getcontext().prec = 4
        list1=os.listdir(dir0);
        cmd1="dcmdump "+list1[0];
        xspace=Decimal(float(str(subprocess.check_output(cmd1+" | grep PixelSpacing", shell=True)).split("\\\\")[0].split("[")[-1]))
        yspace=str(float(str(subprocess.check_output(cmd1+" | grep SliceThickness", shell=True)).split("]")[0].split("[")[-1])-0.01)  #lara -->+0.2)
        zspace=Decimal(float(str(subprocess.check_output(cmd1+" | grep PixelSpacing", shell=True)).split("\\\\")[0].split("[")[-1]))
        Rows=Decimal(str(subprocess.check_output(cmd1+" | grep Rows", shell=True)).split(' ')[2])


        dir2=sequence2[0];
#        for c in dir2:
#                if c=="\\" :
#                    dir2=dir2.replace(c,'')
        dir2=dir2.replace('\\','')
        os.chdir(dir2)
        getcontext().prec = 4
        list1=os.listdir(dir2);
        cmd1="dcmdump "+list1[0];
        xspace2=Decimal(float(str(subprocess.check_output(cmd1+" | grep PixelSpacing", shell=True)).split("\\\\")[0].split("[")[-1]))
        Rows2=Decimal(str(subprocess.check_output(cmd1+" | grep Rows", shell=True)).split(' ')[2])

        pixelSeq2=np.round((float(Rows2)*float(xspace2))/float(xspace));
#        if pixelSeq2!=Rows:
#            xspace=(float(Rows2)*float(xspace2))/float(Rows);      
#            zspace=(float(Rows2)*float(xspace2))/float(Rows);                      

                            
        #os.listdir(patPath+"DCE_MRI_NRRD_scaled/"+dir0.split("/")[-1])[0]
        for file in sequence2:
            os.system("rm -r "+patPath+"DCE_MRI_NRRD_resampled/"+file.split("/")[-1])
            os.system("mkdir "+patPath+"DCE_MRI_NRRD_resampled/"+file.split("/")[-1])
            file2=file
#            for c in file2:
#                if c=="\\" :
#                    file2=file2.replace(c,'')
            file2=file2.replace('\\','')
            f = os.listdir(patPath+"DCE_MRI_NRRD/"+file2.split("/")[-1])[0]
            f2=f
#            if "(" in f2:
#               f2=f2.replace('(','\(')
#            if ")" in filename:
#               f2=f2.replace(')','\)')   
            f2=f2.replace('(','\(')
            f2=f2.replace(')','\)') 
            
            print(xspace)
            print(yspace)
            print(zspace)
#            xspace=1.24
#            yspace=3.2
#            zspace=1.24
       #   os.system("crlResampleToIsotropic -x "+str(xspace)+" -y "+str(yspace)+" -z "+str(zspace)+" "+patPath+"DCE_MRI_NRRD/"+file.split("/")[-1]+"/"+f2+" linear "+patPath+"DCE_MRI_NRRD_resampled/"+file.split("/")[-1]+"/"+f2)
            #os.system("crlResampleToIsotropic -x "+str(xspace-Decimal(0.001))+" -y "+str(yspace)+" -z "+str(zspace-Decimal(0.001))+" "+patPath+"DCE_MRI_NRRD/"+file.split("/")[-1]+"/"+f2+" linear "+patPath+"DCE_MRI_NRRD_resampled/"+file.split("/")[-1]+"/"+f2)      
            #if pixelSeq2!=Rows:
            if xspace!=xspace2:  
                #os.system("crlResampleToIsotropic -x "+str(xspace)+" -y "+str(yspace)+" -z "+str(zspace)+" "+patPath+"DCE_MRI_NRRD/"+file.split("/")[-1]+"/"+f2+" linear "+patPath+"DCE_MRI_NRRD_resampled/"+file.split("/")[-1]+"/"+f2) 
                os.system("crlResampleToIsotropic -x "+str(xspace)+" -y "+str(yspace)+" -z "+str(zspace)+" "+patPath+"DCE_MRI_NRRD/"+file.split("/")[-1]+"/"+f2+" linear "+patPath+"DCE_MRI_NRRD_resampled/"+file.split("/")[-1]+"/"+f2)  
            else:
                os.system("cp "+patPath+"DCE_MRI_NRRD/"+file.split("/")[-1]+"/"+f2+" "+patPath+"DCE_MRI_NRRD_resampled/"+file.split("/")[-1]+"/"+f2) 
         
            #sequencePrefixScaled="COR_RADIAL_VIBE_HI_RES_INIT_T0_12"
            #fileNameScaled=($(find "${folderscaled}/${sequencePrefixScaled}/" -type f -name "*.nrrd"))
            if (fixFOV=='True'):
                print('geometry')
                fileNameScaled=glob.glob(patPath+"DCE_MRI_NRRD_scaled/"+dir0.split("/")[-1]+"/*")[0]
                fileNameScaled2=fileNameScaled
#                for c in fileNameScaled2:
#                    if c=="(":
#                        fileNameScaled2=fileNameScaled2.replace(c,'\(')
#                    if c==")":
#                        fileNameScaled2=fileNameScaled2.replace(c,'\)')
                fileNameScaled2=fileNameScaled2.replace('(','\(')
                fileNameScaled2=fileNameScaled2.replace(')','\)')
                os.system("crlResampler2 -g "+fileNameScaled2+" -i "+patPath+"DCE_MRI_NRRD_resampled/"+file.split("/")[-1]+"/"+f2+" -o "+patPath+"DCE_MRI_NRRD_resampled/"+file.split("/")[-1]+"/"+f2)    
            os.system("crlOrientImage "+patPath+"DCE_MRI_NRRD_resampled/"+file.split("/")[-1]+"/"+f2+" "+patPath+"DCE_MRI_NRRD_resampled/"+file.split("/")[-1]+"/"+f2+" coronal")
        print("Successfully completed resampling of the second sequence.")

#        matrixSize=str(subprocess.check_output(cmd1+" | grep AcquisitionMatrix", shell=True)).split("\\\\")[1]
#        pixelSpacing=str(int(float(str(subprocess.check_output(cmd1+" | grep PixelSpacing", shell=True)).split("\\\\")[0].split("[")[-1])*100)/100)
#        numofSlices=len(os.listdir(dir0))

        
        print("\nConverting NRRD to NIFTY...")
        tNdex=-1;
        for f in sequence2:
            tNdex=tNdex+1;
            nrddFileName=os.listdir(patPath+"DCE_MRI_NRRD_resampled/"+sequence2_0[tNdex].split("/")[-1]);
            nrddFileName1=nrddFileName[0];
#            if "(" in nrddFileName1:
#               nrddFileName1=nrddFileName1.replace('(','\(')
#            if ")" in nrddFileName1:
#               nrddFileName1=nrddFileName1.replace(')','\)')
            nrddFileName1=nrddFileName1.replace('(','\(')
            nrddFileName1=nrddFileName1.replace(')','\)')
            nifPath="/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/"+patientName+"_seq2/";
            os.system("mkdir "+nifPath)
            os.system("crlConvertBetweenFileFormats -in "+patPath+"DCE_MRI_NRRD_resampled/"+f.split("/")[-1]+"/"+nrddFileName1+" -out "+nifPath+"reconSCAN_T"+str(tNdex)+".nii")

        print("Successfully converted sequence 2 from NRRD to Nii.")              
        
        
        
        
#        #ConvertNrrdtoDCEDICOM
#       # if scaleFactor!=1:
#        print("\nConverting NRRD files of the first sequence to DICOM...")
#        folderScaled=patPath+"DCE_MRI_NRRD_scaled/"
#        folderdicom=patPath+"DCE_MRI_corrected/"
#        folder=patPath+"DCE_MRI/"
#        easy=0
#        for file in sequence1:
#                print("Process: "+str(int(easy/len(sequence1)*100))+"%")
#                easy+=1
#                file2=file
#                for c in file2:
#                    if c=="\\" :
#                        file2=file2.replace(c,'')
#
#                if os.path.isdir(patPath+"DCE_MRI_corrected/"+file2.split("/")[-1]):
#                    os.system("rm -r "+patPath+"DCE_MRI_corrected/"+file.split("/")[-1])
#                os.system("mkdir "+patPath+"DCE_MRI_corrected/"+file.split("/")[-1])
#                seqId=file2.split("/")[-1].split("_")[-1]
#                dirName2=file2.split("/")[-1]
#                dirName=file.split("/")[-1]
#                #os.system("mkdir "+patPath+"DCE_MRI_corrected/"+dirName)
#                baseName=os.listdir(patPath+"DCE_MRI_NRRD_scaled/"+dirName2)[0] 
#                baseName2=baseName
#                for c in baseName2:
#                    if c=="(":
#                        baseName2=baseName2.replace(c,'\(')
#                    if c==")":
#                        baseName2=baseName2.replace(c,'\)')
#                dataTemplate=patPath+"DCE_MRI/"+dirName2
#                os.chdir(dataTemplate)
#                list1=os.listdir(dataTemplate);
#                cmd1="dcmdump "+list1[0];
#                study_id=str(subprocess.check_output(cmd1+" | grep StudyID", shell=True)).split("[")[1].split("]")[0]
#                study_id_num=str(subprocess.check_output(cmd1+" | grep StudyInstanceUID",shell=True)).split("]")[0].split("[")[1]
#
#                os.system(rescaleRealToShort+" -s 1 -i "+folderScaled+dirName+"/"+baseName2+" -o "+folderScaled+dirName+"/"+dirName+"_rescaled.nrrd")
#                dataTemplate2=patPath+"DCE_MRI/"+dirName
#                os.system("crlExportDICOM -d "+folderdicom+dirName+" "+folderScaled+dirName+"/"+dirName+"_rescaled.nrrd "+folderdicom+dirName+"/crl_ "+dataTemplate2)
#                os.chdir(folderdicom+dirName2)
#                os.system("cp *0000.dcm seq_template.dcm")
#                os.system("dcmodify -dc -gse seq_template.dcm ")    
#                os.system("rm -rf seq_template.dcm.bak")
#                seq_template_name="seq_template.dcm"    
#                seq_template_name_num=str(subprocess.check_output("dcmdump "+seq_template_name+" | grep SeriesInstanceUID", shell=True)).split("]")[0].split("[")[1]
#                os.system("rm -rf  seq_template.dcm")
#                sliceIdx=0
#                for fileName in glob.glob("*.dcm"):
#                    sliceIdx=int(fileName[-8:-4])
#                    sliceIdxFlip=numofSlices-sliceIdx+1
#                    seriesTime=str(subprocess.check_output("dcmdump "+dataTemplate2+"/"+os.listdir(dataTemplate)[0]+" | grep SeriesTime", shell=True)).split("[")[1].split("]")[0]
#                    acquisitionTime=str(subprocess.check_output("dcmdump "+dataTemplate2+"/"+os.listdir(dataTemplate)[0]+" | grep AcquisitionTime", shell=True)).split("[")[1].split("]")[0]
#                    os.system("dcmodify -i \"(0020,0013)="+str(sliceIdxFlip)+"\" -i \"(0020,0011)="+str(seqId)+"\" -i \"(0020,000e)="+seq_template_name_num+"\" -i \"(0020,0010)="+study_id+"\" -i \"(0020,000d)="+study_id_num+"\" "+fileName)
#                    os.system("dcmodify -i \"(0008,0070)=SIEMENS\" "+fileName)
#                    os.system("dcmodify -i \"(0008,0031)="+seriesTime+"\" -i \"(0008,0032)="+acquisitionTime+"\" "+fileName)
#                os.system("rm -rf *.bak")
#
#        #ChangeTagsConvertNrrdtoDCEDICOM
#
#        print("\nConverting second sequence from NRRD to DICOM...")
#        #os.chdir(sequence1[0])
#        folderresampled=patPath+"DCE_MRI_NRRD_resampled/"
#        folderdicom=patPath+"DCE_MRI_corrected/"
#        folder=patPath+"DCE_MRI/"
#        easy=0
#        for file in sequence2:
#            print("Process: "+str(int(easy/len(sequence2)*100))+"%")
#            easy+=1
#            file2=file
#            for c in file2:
#                if c=="\\" :
#                    file2=file2.replace(c,'')
#            if os.path.isdir(patPath+"DCE_MRI_corrected/"+file2.split("/")[-1]):
#                os.system("rm -r "+patPath+"DCE_MRI_corrected/"+file.split("/")[-1])
#            os.system("mkdir "+patPath+"DCE_MRI_corrected/"+file.split("/")[-1])
#            dirName2=file2.split("/")[-1]
#            dirName=file.split("/")[-1]
#            seqId=file2.split("/")[-1].split("_")[-1]
#
#            baseName=os.listdir(patPath+"DCE_MRI_NRRD_resampled/"+dirName2)[0]
#            baseName2=baseName
#            for c in baseName2:
#                if c=="(":
#                    baseName2=baseName2.replace(c,'\(')
#                if c==")":
#                    baseName2=baseName2.replace(c,'\)')
#            dataTemplate=folder+dirName2
#            os.chdir(dataTemplate)
#            list1=os.listdir(dataTemplate);
#            cmd1="dcmdump "+list1[0];
#            study_id=str(subprocess.check_output(cmd1+" | grep StudyID", shell=True)).split("[")[1].split("]")[0]
#            study_id_num=str(subprocess.check_output(cmd1+" | grep StudyInstanceUID",shell=True)).split("]")[0].split("[")[1]
#            os.system(rescaleRealToShort+" -s 1 -i "+folderresampled+dirName+"/"+baseName2+" -o "+folderresampled+dirName+"/"+dirName+"_rescaled.nrrd")
#            dataTemplate2=patPath+"DCE_MRI/"+dirName
#            os.system("crlExportDICOM -d "+folderdicom+dirName+" "+folderresampled+dirName+"/"+dirName+"_rescaled.nrrd "+folderdicom+dirName+"/crl_ "+dataTemplate2)
#            os.chdir(folderdicom+dirName2)
#            os.system("cp *0000.dcm seq_template.dcm")
#            os.system("dcmodify -dc -gse seq_template.dcm ")    
#            os.system("rm -rf seq_template.dcm.bak")
#            seq_template_name="seq_template.dcm"    
#            seq_template_name_num=str(subprocess.check_output("dcmdump "+seq_template_name+" | grep SeriesInstanceUID", shell=True)).split("]")[0].split("[")[1]
#            os.system("rm -rf  seq_template.dcm")
#            sliceIdx=0
#            for fileName in glob.glob("*.dcm"):
#                sliceIdx=int(fileName[-8:-4])
#                sliceIdxFlip=numofSlices-sliceIdx+1
#                seriesTime=str(subprocess.check_output("dcmdump "+dataTemplate2+"/"+os.listdir(dataTemplate)[0]+" | grep SeriesTime", shell=True)).split("[")[1].split("]")[0]
#                acquisitionTime=str(subprocess.check_output("dcmdump "+dataTemplate2+"/"+os.listdir(dataTemplate)[0]+" | grep AcquisitionTime", shell=True)).split("[")[1].split("]")[0]         
#                os.system("dcmodify -i \"(0020,0013)="+str(sliceIdxFlip)+"\" -i \"(0020,0011)="+str(seqId)+"\" -i \"(0020,000e)="+seq_template_name_num+"\" -i \"(0020,0010)="+study_id+"\" -i \"(0020,000d)="+study_id_num+"\" "+fileName)
#                os.system("dcmodify -i \"(0018,1310)=0\\"+matrixSize+"\\"+matrixSize+"\\0"+"\" -i \"(0028,0010)="+matrixSize+"\" -i \"(0028,0011)="+matrixSize+"\" -i \"(0028,0030)="+pixelSpacing+"\\"+pixelSpacing+"\" "+fileName)
#                os.system("dcmodify -i \"(0008,0070)=SIEMENS\" "+fileName)
#                os.system("dcmodify -i \"(0008,0031)="+seriesTime+"\" -i \"(0008,0032)="+acquisitionTime+"\" "+fileName)
#            os.system("rm -rf *.bak")
#        print("Successfully converted second sequence from NRRD to DICOM.")
#
#
#        for direc in glob.glob(patPath+"*"+"DCE_MRI_corrected/*"):
#            os.chdir(direc)
#            os.system("rm *.bak")
    else:
        print(numSeq)
        sys.exit("The number of sequences was entered wrong. The acceptable values are 1 and 2.")
    os.system("chmod 774 -R "+patPath[:-(len(patPath.split("/")[-2])+1)])
    print("\n\tThe data is ready to use.")
    
    
    
def dicomHeaderInfoExtractor(patientName,patPath, seqPrefix1, seqPrefix2,numSeq):

#############################################  

    if not patPath.endswith("/"):
        patPath=patPath+"/"
    files=glob.glob(patPath+"*"+seqPrefix1+"*")
    files=sorted(files, key=lambda name: int(name[len(patPath):].split('_')[0]));
    if len(files)==0:
        sys.exit("No directories with the given sequence prefix have been found.")

    # read sequences 
    sequence1=[];sequence1Pcorrected=[];
    for filename in files:
        sequence1.append(filename)
#        if "(" in filename:
#           filename=filename.replace('(','\(')
#        if ")" in filename:
#           filename=filename.replace(')','\)')      
        filename=filename.replace('(','\(')
        filename=filename.replace(')','\)')      
        sequence1Pcorrected.append(filename)

    tNdex=-1;
    sec1=[];
    for f in sequence1:
        tNdex=tNdex+1;
        dicomFileName=os.listdir(patPath+"/"+sequence1[tNdex].split("/")[-1]);
    #         print(dicomFileName1);
        cmd1="dcmdump "+sequence1Pcorrected[tNdex]+'/'+dicomFileName[0];
        sn=str(subprocess.check_output(cmd1+" | grep StationName", shell=True)).split("[")[1].split("]")[0]
        ID=str(subprocess.check_output(cmd1+" | grep PatientID", shell=True)).split("[")[1].split("]")[0]
        age=str(subprocess.check_output(cmd1+" | grep PatientAge", shell=True)).split("[")[1].split("]")[0]
        mm=str(subprocess.check_output(cmd1+" | grep ManufacturerModelName", shell=True)).split("[")[1].split("]")[0]

        time=str(subprocess.check_output(cmd1+" | grep AcquisitionTime", shell=True)).split("[")[1].split("]")[0]
        sec=int(time.split(".")[0][0:2])*3600+int(time.split(".")[0][2:4])*60+int(time.split(".")[0][4:6])+(float(time)-int(float(time)))
        sec1.append(sec);
    
    seq1VolNum=len(sequence1);  
    seq2VolNum=0;  
    seq1EndTime=sec1[-1];   
    seq2StartTime=0;
    timeRes0=np.diff(np.array(sec1),n=1);
    if np.size(timeRes0)==0:
        timeRes=[];
    else:
        timeRes=timeRes0[0];
#StationName
    ############################################# if num of sequence is 2   
    if numSeq==2:
        if not patPath.endswith("/"):
            patPath=patPath+"/"
        files=glob.glob(patPath+"*"+seqPrefix2+"*")
        files=sorted(files, key=lambda name: int(name[len(patPath):].split('_')[0]));
        if len(files)==0:
            sys.exit("No directories with the given sequence prefix have been found.")

        # read sequences 
        sequence1=[];sequence1Pcorrected=[];
        for filename in files:
            sequence1.append(filename)
#            if "(" in filename:
#               filename=filename.replace('(','\(')
#            if ")" in filename:
#               filename=filename.replace(')','\)')    
            filename=filename.replace('(','\(')
            filename=filename.replace(')','\)') 
            sequence1Pcorrected.append(filename)

        tNdex=-1;
        sec1=[];
        for f in sequence1:
            tNdex=tNdex+1;
            dicomFileName=os.listdir(patPath+"/"+sequence1[tNdex].split("/")[-1]);
    #         print(dicomFileName1);
            cmd1="dcmdump "+sequence1Pcorrected[tNdex]+'/'+dicomFileName[0];
            time=str(subprocess.check_output(cmd1+" | grep AcquisitionTime", shell=True)).split("[")[1].split("]")[0]
            sec=int(time.split(".")[0][0:2])*3600+int(time.split(".")[0][2:4])*60+int(time.split(".")[0][4:6])+(float(time)-int(float(time)))
            sec1.append(sec);

        seq2VolNum=len(sequence1);   
        seq2StartTime=sec1[0];   
        timeRes0=np.diff(np.array(sec1),n=1);
        timeRes2=timeRes0[0];

        timeRes=[timeRes,timeRes2];


    dicomHeaderInfo= {'timeRes':timeRes, 'seq1EndTime':seq1EndTime, 'seq2StartTime':seq2StartTime, 'StationName':sn,\
                      'ManufacturerModelName':mm, 'PatientID':ID,'PatientAge':age,'seq1VolNum':seq1VolNum,'seq2VolNum':seq2VolNum};
    return dicomHeaderInfo    


