%%
%
%       This script processes the subjects in the MRU list.
%
%       PROCESSING STEPS:
%           - Load k-space data
%                   - loads the kspace data (k22n)
%                   - computes the FFT in the stack direction (k3n)
%                   - computes the position of k-space samples  (k_samples)
%                   - computes the compensation function for the k-space samples (dcf)
%           - downsamples the data in k-space by eliminating the samples further from the center
%                   - selects (fraction_of_points) of the overall data
%                   - resamples all relevant files and recomputes the new sizes
%           - compute coil profiles
%                   - computes the coil profiles of all samples
%           - compute FIDs
%                   - most importantly computes the correlation based metric
%           - computes NUFFT (no-regularization) reconstructions
%                   - joins all the results in single 4D file too
%
%
%
%
clc;close all; clear all;

addpath(genpath('./ext/')); 
%removed - obsolete - run 'addpaths.m'; %sergev - instructed by Jaume to add this as the start of recon 
parpool('local',4);

%% initialize run parameters
run_parameters = struct('script','process_all_subjects.m');

% output directory 
root_data_folder = './output/'; % save to the same directory where the scripts is located 
    % > CEMRE - Sila asked for output to be written here but I dont have write access - /fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsGRASP/method2/     


% directory with SCANNER reconstructions (aka raw data was reconstructed ON the scanner with Siemens software and copied to the following directory) 
reconResultsSCANdir = '/fileserver/external/rdynaspro4/abd/GraspRecons/reconResultsSCAN/'; % sergev - corrected the reconResultsSCAN path. old path - '/common/abd/GraspRecons/reconResultsSCAN/'; 

% FOR FLYWHEEL ONLY - setting a temporary directory to hold example data 
reconResultsSCANdir = './input/reconResultsSCAN/'; % sergev - corrected the reconResultsSCAN path. old path - '/common/abd/GraspRecons/reconResultsSCAN/'; 


%% PARAMS
patientName = {'Name_Surname'};   % copy and paste correct value from `../input/subjectList_MRU.csv`. Remember to never store PHI data on github.   
                % 'Name_Surname2',...  % possible to add more subjects to process at the same time           
                % 'Name_Surname3'};  
run_parameters.patientName = patientName;
run_parameters.subjectListCSV = '../input/subjectList_MRU.csv'; %sergev - added local .csv with `END` written in each cell - else matlab cannot find line end during reading process

run_parameters.spokes_per_vol = 34;        % number of stokes per volume to use (minimal is 34)
run_parameters.fraction_of_points = 1;   % number of points per line to use in the reconstuction
run_parameters.lambda_GRASP = 0.0125;      % regularization prameter for GRASP reconstruction
run_parameters.FIDcorr_th = 0.01;          % threshold to consider corrupted line in correlation FIDS
run_parameters.percentage_th = 0.5;        % percentage of spokes in a volume to consider it corrupted
run_parameters.contrastTime = 20;          % initial time without contrast. (in secs and assuming ~100ms per spoke) 
run_parameters.coil_err_th = 0.75;         % threshold to select bad coils
contrastSpokes = 1;%contrastTime/0.1;

%% retrieve patient info
[  files, timePoints, patientName , MRN, referenceNii, badCoils] = parseCSVforSubjectData( patientName, run_parameters.subjectListCSV );

%% start report
textname = './output/report.txt'; 
fo = fopen(textname,'w+'); %sergev - writing a report in local directory 

%% iterate through al subjects
for ss = 1:numel(files)

    reportLine = sprintf('SUBJECT: %d, %s. (%d of %d)\n', MRN{ss}, patientName{ss},  ss, numel(files)); 
    fprintf(reportLine);
    fprintf(fo,reportLine);

    %% if no raw files found ignore
    if numel(files{ss}) > 0

        %% create folders
        subject_processed_folder = sprintf('%s/DCE/Processed/Subj%d/',root_data_folder,MRN{ss});
        mkdir(subject_processed_folder);
        subject_processed_folder = sprintf('%s/%s/',subject_processed_folder, patientName{ss});
        mkdir(subject_processed_folder);
        subject_processed_folder = sprintf('%s/common_processed/', subject_processed_folder );
        mkdir(subject_processed_folder);
        

        %% load data
        [  k22n, k3n, k_samples, dcf, totalSpokesInSeq, SliceOversampling, coilList] = load_kspace_data( files{ss}, timePoints{ss}, run_parameters.spokes_per_vol ); % sergev - added 'magnetom vida' option for loading raw data
        [N,NCha, NSpk,NSli] = size(k3n);
        

        %% extract FIDs
        FID_folder = sprintf('%s/FID_folder_fast/',subject_processed_folder);
        FID_filename = sprintf('%s/subj%d_FID_signal.mat',FID_folder, MRN{ss});
        FIDcorr_filename = sprintf('%s/subj%d_FID_corr_avrg_signal.mat',FID_folder, MRN{ss});
        FIDall_filename = sprintf('%s/subj%d_FID_all_avrg_signal.mat',FID_folder, MRN{ss});
        if ~exist(FID_filename) | ~exist(FIDcorr_filename)  | ~exist(FIDall_filename)

            reportLine = sprintf('Running FIDs...\n');
            fprintf(reportLine);
            fprintf(fo,reportLine);

            mkdir(FID_folder);

            % extract FIDs
            [FID_signal,FID_corr, FID_meanNormChange, FID_maxNormChange, FID_medNormChange] = FID_extraction(k22n);
            save(FID_filename,'FID_signal');
            save(FIDcorr_filename,'FID_corr');
            save(FIDall_filename,'FID_signal','FID_corr','FID_meanNormChange', 'FID_maxNormChange', 'FID_medNormChange');

        else
            reportLine = sprintf('%s already exists\n',FID_filename);
            fprintf(reportLine);
            fprintf(fo,reportLine);
            load(FIDall_filename);
        end

        clear 'k22n';


        %% select bad coils
        coil_folder = sprintf('%s/coil_profile_N%d/',subject_processed_folder, N);
        coil_images = sprintf('%s/subj%d_coilprofile_N%d_rec_images.mat',coil_folder, MRN{ss}, N);
        if ~exist(coil_images)

            mkdir(coil_folder);

            % prepare parallel loop
            k3cell = cell(1,NCha);
            ksampCell = cell(1,NCha);
            dcfCell = cell(1,NCha);
            for ch = 1:NCha
                k3cell{ch} = k3n(:,ch,1:34*6,:);
                ksampCell{ch} = k_samples(:,1:34*6,:);
                dcfCell{ch} = dcf(:,1:34*6,:);
            end

            coil_rec_image = cell(1,NCha);
            parfor ch = 1:NCha
                [coil_rec_image{ch}] = image_reconstruction_noregu( k3cell{ch}, ksampCell{ch}, dcfCell{ch}, ones(N,N,1,NSli) , run_parameters.spokes_per_vol*6);
            end
            tempCoil = zeros( N,N,NSli, NCha );
            for ch = 1:NCha
                tempCoil(:,:,:,ch) = coil_rec_image{ch}{1};
            end
            coil_rec_image = tempCoil;
            save(coil_images,'coil_rec_image','-v7.3');

            clear tempCoil k3cell ksampCell dcfCell            

        else
            fprintf('%s already exists\n',coil_images);
            load(coil_images);
        end

        % select bad coils
        [goodCoils, covErr] = selectBadCoils( coil_rec_image, [10 10], [100,100], run_parameters.coil_err_th );

        clear('coil_rec_image');    % done with the coil images
        
        run_parameters.goodCoils = goodCoils;

        %% remove bad stokes and coils
        lowres_k_samples = [  floor(-N*run_parameters.fraction_of_points/2+1) : ceil(N*run_parameters.fraction_of_points/2)  ] + floor(N/2);

        goodCoilsIX = find(goodCoils == 1);

        k3n = k3n(lowres_k_samples,goodCoilsIX,:,:);
        k_samples = k_samples(lowres_k_samples,:,:) * 1/run_parameters.fraction_of_points;  % the multiplication compensates a weird scaling
        dcf = dcf(lowres_k_samples,:,:);

        [N,NCha, NSpk,NSli] = size(k3n);


        %% compute coil profiles
        coil_folder = sprintf('%s/coil_profile_N%d_autoBadCoils/',subject_processed_folder, N);
        coil_filename = sprintf('%s/subj%d_coilprofile_N%d.mat',coil_folder, MRN{ss}, N);
        if ~exist(coil_filename)

           reportLine = sprintf('Running Coil Profiles...\n');
           fprintf(reportLine);   
           mkdir(coil_folder);
           
           k3cell = cell(1,NSli);
           for sl = 1:NSli
               k3cell{sl} = k3n(:,:,1:run_parameters.spokes_per_vol*6,sl);
           end
           
           coilprofile = cell(1,NSli);
           parfor sl = 1:NSli
               coilprofile{sl} = calc_iterativeradial_multicoil_onechunk( 100*k3cell{sl} );
           end
           coilprofile = cell2mat( reshape( coilprofile, [1,1,1,NSli] ) );
           save(coil_filename,'coilprofile', '-v7.3');
           
           clear k3cell

        else
            fprintf('%s already exists\n',coil_filename);
            load(coil_filename);
        end

        

        %% compute GRAST and NONREGU reconstructions
        GRASP_folder = sprintf('%s/GRASP_reconstruction_autoBadCoils/',subject_processed_folder);
        GRASP_folder_OLD = sprintf('%s/../',subject_processed_folder);
        GRASP_filename = sprintf('subj%d_recon_GRASP_N%d_S%d_lambda%0.6f_vol4D',MRN{ss},N, run_parameters.spokes_per_vol, run_parameters.lambda_GRASP);
        GRASP_filename_OLD = sprintf('recon_GRASP_4D');
        
        
        GRASP_metadata_filename = sprintf('subj%d_recon_GRASP_N%d_S%d_lambda%0.6f_metadata.mat',MRN{ss},N, run_parameters.spokes_per_vol);               
        if ~exist([GRASP_folder,GRASP_filename,'.mat'])

            if numel(referenceNii{ss}) == 0
                if length(timePoints{ss})>1
                    referenceNii{ss} = strcat(reconResultsSCANdir,patientName{ss},'_seq1/reconSCAN_T0.nii');
                else
                    referenceNii{ss} = strcat(reconResultsSCANdir,patientName{ss},'/reconSCAN_T0.nii');
                end
            end

            reportLine = sprintf('Running NUFFT reconstruction with N%d points and %d stokes...\n', N, run_parameters.spokes_per_vol);
            fprintf(reportLine);

            % make dir
            mkdir(GRASP_folder);
            [rec_image_GRASP] = image_reconstruction_wrapper_GRASP( k3n,  k_samples, dcf, coilprofile, run_parameters.spokes_per_vol, run_parameters.lambda_GRASP, contrastSpokes);

            % writeGraspReconMat2nii_jcf(recon_GRASP, patientName{ss}, GRASP_folder, GRASP_filename, timePoints{ss}, SliceOversampling);
            save( sprintf( '%s/%s.mat', GRASP_folder, GRASP_filename) , 'rec_image_GRASP', 'run_parameters', '-v7.3');
%             save( sprintf( '%s/%s.mat', GRASP_folder, GRASP_metadata_filename), 'git_HEAD_version', 'goodCoilsIX', 'lambda_GRASP','spokes_per_vol', 'contrastTime','fraction_of_points' );
            
            %% save as nii
            numImgs = size(rec_image_GRASP,3);
            cellImg = cell(1,numImgs);
            for tt = 1:numImgs
                tempImg = zeros( N,N,NSli );
                tempImg = squeeze(rec_image_GRASP(:,:,tt,:));
                cellImg{tt} = tempImg;
            end
            clear tempImg;
            writeGraspReconMat2nii_jcf(cellImg, patientName{ss}, GRASP_folder, GRASP_filename, SliceOversampling, referenceNii{ss}, true);
            clear cellImg;

        else
            reportLine = sprintf('%s already exists\n',GRASP_filename);
            fprintf(reportLine);
        end

             clear('k22n', 'k3n', 'coilprofile', 'k_samples', 'dcf', 'totalSpokesInSeq');

    else
        reportLine = sprintf('FAILED SUBJECT: %d, %s. (%d of %d)\n', MRN{ss}, patientName{ss},  ss, numel(files));
        fprintf(reportLine);
        fprintf(fo,reportLine);
    end;

end;

fclose(fo);