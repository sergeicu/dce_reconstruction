############ RUNNING EXAMPLE RECONSTRUCTION ############

# Step 1: open matlab and cd into `reconstruction_matlab/` directory 
cd reconstruction_matlab/

# Step 2: 
# open `process_all_KidneySubjects.m` and substitute the correct `patientName` (line46)  with information from `input/subjectList_MRU.csv`. Remember to never store PHI data on github.


# Step 3: run example 
process_all_Kidney_subjects.m

# NOTES 
# it is recommended that you run reconstruction on a server with lots of memory... e.g. no less than 100Gb RAM 
# so instead of running it on your local machine you can ssh into one of the following CRL machines and run the reconstruction. e.g 
ssh -X boreas #or ssh ganymede / auster / zephyr / io 
# make sure that you run with '-X' command so that you can start matlab GUI within each machine. 




############ REQUIRED CHANGES FOR PROCESSING NEW RECONSTRUCTIONS ############

# 1. Search for `FOR FLYWHEEL ONLY` text inside `process_all_KidneySubjects.m`. The note should be pretty self explanatory. i.e. just comment this line out - and it would default to hospital directory. 

# 2. Search for `CEMRE` text inside `process_all_KidneySubjects.m`. The note should be pretty self-explanatory. In short - ideally, we should store all new reconstructions to the correct output folder. However, I could not write to it as I don't have permission. You need to ask Sila to give you write permission to that folder. 

# 3. If you want to run reconstruction on a NEW subject you need to do the following: 

# 3.1. Change `patientName` in `process_all_KidneySubjects.m`

# 3.2. Add patient name, MRN and raw data location to `../input/subjectList_MRU.csv`

# 3.3. Check that `reconResultsSCANdir` points to the correct folder. This folder must have the new patient's scanner reconstruction (exported from the siemens scanner) inside it. 





############ To change the output file naming convention to the old format ############


# If you want to name the files according to the old format -> open the matlab function and substitute the following: 
#    inside these lines: 
#    `writeGraspReconMat2nii_jcf(cellImg, patientName{ss}, GRASP_folder, GRASP_filename, SliceOversampling, referenceNii{ss}, true);`
#    and 
#    `save( sprintf( '%s/%s.mat', GRASP_folder, GRASP_filename) , 'rec_image_GRASP', 'run_parameters', '-v7.3');`
#    Change: 
#    GRASP_folder --> GRASP_folder_OLD 
#    GRASP_filename --> GRASP_filename_OLD 

#    also comment out this line 
#    `mkdir(GRASP_folder);` -> LINE 231


#     FINALLY, EXTREMELY IMPORTANT: 
#    you must COMMENT OUT the following `if` statement completely (otherwise the script will never proceed to compute the images since the GRASP_folder already exists and the if statement would skip all of the code inside it (THIS ONLY HAPPENS IF YOU WILL USE `GRASP_folder_OLD`) 
#    `if ~exist([GRASP_folder,GRASP_filename,'.mat'])




#    FINAL COMMENT: 
#    warning: i haven't really tested reconstruction with the OLD file format. This is because I would need to re-run an entire reconstruction pipeline again (and it would take 24 hours of time). So you may need to debug more. 

#   TO COMPARE THE OLD FORMAT TO NEW FORMAT CHOOSE THE FOLLOWING FOLDERS:  
# 
#    OLD FORMAT
#    `/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsGRASP/method2/`   # find a recent subject here 
#    
#    NEW FORMAT 
#    `/fileserver/abd/serge/dce_data/reconstructed_matlab/`  # two of the subjects here have the correct format of reconstruction 