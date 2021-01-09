

function [] = loadAndSaveKspaceData( patientName, subjectListCSV, spokes_per_vol, outputFile, subject_processed_folder )

    % params
    coilEstimationSpokes = 1:spokes_per_vol*6;
    coil_err_th = 0.75;         % threshold to select bad coils

    % parse information from the CSV file
    disp('gato')
    [  files, timePoints, patientName , MRN, referenceNii, badCoils] = parseCSVforSubjectData( patientName, subjectListCSV );

files

    % load data
    disp('gato2')
    files{1}
    [  k22n, k3n, k_samples, dcf, totalSpokesInSeq, SliceOversampling, coilList] = load_kspace_data( files{1}, timePoints{1}, spokes_per_vol );

    [N,NCha, NSpk,NSli] = size(k3n);

    % run FID for outlier detection
    FID_folder = sprintf('%s/FID_folder_fast/',subject_processed_folder);
    FID_filename = sprintf('%s/subj%d_FID_signal.mat',FID_folder, MRN{1});
    FIDcorr_filename = sprintf('%s/subj%d_FID_corr_avrg_signal.mat',FID_folder, MRN{1})
    FIDall_filename = sprintf('%s/subj%d_FID_all_avrg_signal.mat',FID_folder, MRN{1})
    if ~exist(FID_filename) | ~exist(FIDcorr_filename)  | ~exist(FIDall_filename)

        reportLine = sprintf('Running FIDs...\n');
        fprintf(reportLine);

        mkdir(FID_folder);

        % extract FIDs
        [FID_signal,FID_corr, FID_meanNormChange, FID_maxNormChange, FID_medNormChange] = FID_extraction(k22n);
        save(FID_filename,'FID_signal');
        save(FIDcorr_filename,'FID_corr');
        save(FIDall_filename,'FID_signal','FID_corr','FID_meanNormChange', 'FID_maxNormChange', 'FID_medNormChange');

    end

    % run FID to estimate the respiratory phase
    FID_respPhase_filename = sprintf('%s/subj%d_FID_respPhaseSignal.mat',FID_folder, MRN{1});
    if ~exist(FID_respPhase_filename) 

        reportLine = sprintf('Running estimation of respiratory phase...\n');
        fprintf(reportLine);

        [Res_Signal] = FID_respPhase_estimation( k22n );
        save(FID_respPhase_filename,'Res_Signal');

    end
    
    % run coilprofiles
    coil_folder = sprintf('%s/coil_profile_N%d/',subject_processed_folder, 448);
    coil_filename = sprintf('%s/subj%d_coilprofile_N%d.mat',coil_folder, MRN{1},  448);
    coil_all_filename = sprintf('%s/subj%d_coilprofile_allcoils_N%d.mat',coil_folder, MRN{1},  448);
    badCoils_filename = sprintf('%s/subj%d_badCoils.mat',coil_folder, MRN{1});
    if ~exist(coil_filename)

        reportLine = sprintf('Running Coil Profiles...\n');
        fprintf(reportLine);   
        mkdir(coil_folder);

        % compute all coils
        k3cell = cell(1,NSli);
        for sl = 1:NSli
            k3cell{sl} = k3n(:,:,coilEstimationSpokes,sl);
        end
        coilprofile = cell(1,NSli);
        parfor sl = 1:NSli
            coilprofile{sl} = calc_iterativeradial_multicoil_onechunk( 100*k3cell{sl} );
        end
        coilprofile = cell2mat( reshape( coilprofile, [1,1,1,NSli] ) );
        save(coil_all_filename,'coilprofile', '-v7.3');
        clear coilprofile


        % prepare data
        k3cell = cell(1,NCha);
        ksampCell = cell(1,NCha);
        dcfCell = cell(1,NCha);
        for ch = 1:NCha
            k3cell{ch} = k3n(:,ch,coilEstimationSpokes,:);
            ksampCell{ch} = k_samples(:,coilEstimationSpokes,:);
            dcfCell{ch} = dcf(:,coilEstimationSpokes,:);
        end
     
        % determine bad coils
        coil_rec_image = cell(1,NCha);
        parfor ch = 1:NCha
            [coil_rec_image{ch}] = image_reconstruction_noregu( k3cell{ch}, ksampCell{ch}, dcfCell{ch}, ones(N,N,1,NSli) , spokes_per_vol);
        end
        
        tempCoil = zeros( N,N,NSli, NCha );
        for ch = 1:NCha
            tempCoil(:,:,:,ch) = coil_rec_image{ch}{1};
        end

        [goodCoils, covErr] = selectBadCoils( tempCoil, [10 10], [100,100], coil_err_th );
        goodCoilsIX = find(goodCoils == 1);
        save(badCoils_filename, 'goodCoilsIX')

        clear tempCoil k3cell ksampCell dcfCell  coil_rec_image    

        % eliminate bad coils from data
        k3n = k3n(:,goodCoilsIX,:,:);

        % compute coils
        k3cell = cell(1,NSli);
        for sl = 1:NSli
            k3cell{sl} = k3n(:,:,coilEstimationSpokes,sl);
        end
        
%         parpool('local',4);
        coilprofile = cell(1,NSli);
        parfor sl = 1:NSli
            coilprofile{sl} = calc_iterativeradial_multicoil_onechunk( 100*k3cell{sl} );
        end
        coilprofile = cell2mat( reshape( coilprofile, [1,1,1,NSli] ) );
        save(coil_filename,'coilprofile', '-v7.3');
        
        clear k3cell

    else
        fprintf('%s already exists\n',coil_filename);

        if exist(badCoils_filename)
            load(badCoils_filename);
            % eliminate bad coils from data
            k3n = k3n(:,goodCoilsIX,:,:);
        end
    end

    % save k-space data
    save(outputFile, 'k3n', 'k_samples', 'dcf', 'totalSpokesInSeq', 'SliceOversampling', 'coilList', '-v7.3');


end
