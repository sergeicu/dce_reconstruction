PLEASE IGNORE THIS FILE IF YOU ARE RUNNING DEMO DCE RECONSTRUCTION. 

This file contains instructions on data preprocessing for new subjects only. 

## STEP 0: Setup 

Please follow instructions to install new conda environment with python dependencies in `instructions/reconstruction/python_run_example.sh`

Activate your new python environment 
`conda activate py3.7`

You must also make sure that you can access the following file: 
`/fileserver/external/rdynaspro4/abd/MRUcommon/retrieve2.sh` (script for pulling DICOM images from the scanner)


You must also have access to CRKit tools, such as: `crlConvertBetweenFileFormats`. Please refer to CRL-guide on how to set this up. Or ask Serge / Simon for assistance. 

##  STEP 1: Get patient details. 

You will need this information to proceed with preprocessing and reconstruction. 
- FirstName_LastName 
- MRN
- Date of Acquisition
- full path to raw data file (.dat) 

Raw data is typically stored here: 
    - /fileserver/abd/radial
    - /fileserver/external/rdynaspro4/abd/radial_JCF/
    - /fileserver/projects/abd/radial_ralph/
However, it is best that you can Sila for direct full path right away. 

## Step 2: Update .csv and .xls

Update the following local file: 
- `../input/subjectList_MRU.csv` - required for reconstruction and preprocessing   

Update the following CRL files (on the CRL filesystem): 
In:`/fileserver/external/rdynaspro4/abd/MRUcommon/` 
- `subjectList_MRU.xls` 
- `subjectDicomInfo.xls`
.xls update is required for information and historical purposes only (not currently used in preprocessing or reconstruction)


## STEP 3: Check DICOMs

Check if patient DICOM files have already been copied onto the local drive from the scanner. 

Search for existence of the folder with patient name here: 
`/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN`

Note that you may already see patient's name and MRN inside the .xls file as well (which indicates that someone had run the preprocessing pipeline already).

## STEP 4: Edit `singlePatientFullProcessing.py`

Open `singlePatientFullProcessing.py` in a text editor/IDE


Edit FLAGS to enable/disable these steps of the preprocessing pipeline:

- downloadDataEnabled=1 - set to 0 for multiCenterStudy (CHOA or NCH) and run `prepareMulticenterStudy_CHOA.py` or `prepareMulticenterStudy_NCH.py` instead (see separate instructions)
- preprocessingEnabled=1 
- visualizationEnabled=0 - always set to zero (not used currently)
- createMIPs=1
- kidneyMaskEnabled=0 - set to 1 ONLY if you are using the correct environment (see below)
- aortaMaskEnabled=0 - set to 1 ONLY if you are using the correct environment (see below)


[WARNING]
Certain functionalities in the preprocessing pipeline rely on Jaume Coll Font's libraries that are stored *locally* on the CRL's servers. This does NOT apply to reconstruction code (only to preprocessing code). For example, preprocessing flags such as `kidneyMaskEnabled` and `aortaMaskEnabled` require Jaume's libriaries to be enabled. 

Also I did not personally need those enabled to perform basic pre-processing. However, if you want to run these, you will need to export binaries from Jaume's local libraries. 

To do this, please run: 
`source /home/ch199899/miniconda3_bkup_startup.sh` 

And then proceed to enable these flags: 
- kidneyMaskEnabled=1
- aortaMaskEnabled=1


## STEP 4: Run preprocessing 

`cd preprocessing `
`python singlePatientFullProcessing.py `

You will be asked to enter patientName, MRN and date. 
Please enter them in the same format as listed in the .xls / .csv file for other patients. 

Note that the script will stop and output: 
`Execution paused. Delete bad folders`

At this point, you should manually go into `/fileserver/commondataraw/MRUoriginal/<PatientName>` and check that there is only a single subfolder here. Else the script may fail. Note that this check should be automated over time. 

-- 

## RECONSTRUCTION FOR multiCenter data (NCH or CHOA)

We have two large external collaborations for DCE study - these are referred to as NCH and CHOA respectively. 

These studies send us their preprocessed DICOM files as well as raw files (.dat). 
As a result, these files require a different reconstruction pipeline. 


### Step 1: Prepare custom NCH/CHOA scripts
To run multicenter study you need to operate with these scripts: 
- `prepareMultiCenterStudy_CHOA.py`
- `prepareMultiCenterStudy_NCH_serge.py`

Note that you still need to go through each of the steps above (Steps 1-4). E.g. check if the DICOM directory exists (Step 2), etc

Open corresponding .py file (NCH or CHOA) and fill in the correct entries for: 
    - subjectName 
    - mrn 
    - scanDate 
    - folderTag 
    - T 

Note that manually imported DICOM files for multicenter data are typically stored in: 
- `/fileserver/external/rdynaspro4/abd/multiCenterStudy/RAW`

To set `mrn`: 
- set it to some fake number
- e.g. we typically set it to the next number in the sequence (although it is not directly necessary)

To set `scanDate`:
- similar to above - set a fake scan data (or date when the data was imported from NCH / CHOA study)

To set `folderTag`: 
- Check this: 
- ls `/fileserver/external/rdynaspro4/abd/multiCenterStudy/RAW/${subj}/DICOM/`, where ${subj} is the subject named folder
- note the 'folder tag' - i.e. all DCE dicom folders should have a common name convention, that is not present in NON-dce dicom folders (scans) 
- put this as 'folderTag' in the .py

To set `T`:

- count how many DCE volumes there are in the DICOM directory above  
- e.g. 
- `ls /fileserver/external/rdynaspro4/abd/multiCenterStudy/RAW/{subj}/DICOM/ | wc -l`
- put this as `T` parameter in the .py file 

### Step 2: Run custom NCH/CHOA script

Run appropriate script to generate the correct directories in the main folder. 

### Step 3: Manually update the .xls and .csv 

Unfortunately this script does not update the .xls file automatically. Hence you need to download the .xls file manually and update it manually (and upload it back to the server). 

Same for .csv file (stored locally in `../input/subjectList_MRU.csv`)

Furthermore, please read instructions at the bottom of the script. 

### Step 4: 
In `singlePatientFullProcessing.py` set `downloadDataEnabled=0` because the data does NOT exist on the hospital PACS system (i.e. where the original scans are located). Instead, the data exists in the local folder (where Sila had chosen to store the NCH/CHOA DICOM and .dat files). 

### Step 5: 



Next, proceed with matlab/python recon steps. 