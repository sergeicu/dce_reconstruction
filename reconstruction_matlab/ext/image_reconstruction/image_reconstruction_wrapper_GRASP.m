%% 
%
%       This function applies the reconstruction of the MR image from 
%       k-space sampling.
%
%       It uses Multi-Coil Non-Uniform Fast Fourier Transform (MC-NUFFT) for the reconstruction. 
%       Thus, it relies on the nufft_toolbox by Jeff Fessler as well as the single-coil NUFFT 
%       operator from Miki Lustig and the multi-coil operator from Li Feng & Ricardo Otazo.
%
%       This code is a re-wrapping of the implementation that Ali Pour Yazdanpanah did at CRL (2017)
%       which in turn is a re-package of Sila Kurugol's and Marzieh's Haggigis methods.
%
%       INPUTS:
%           - k_data - <N,NCha,NSpk,NSli>double - data sampled in K-space.
%                                                   Fourier transformed along the "cartesian" direction
%                                               Dimensions: - N = number of samples per spoke.
%                                                           - NCha = number of chanels.
%                                                           - NSpk = number of spokes.
%                                                           - NSli = number of slices.
%           - k_samples - <baseresolution,1>double -    k space coordinates ascomplex values (kx + i ky)
%           - dcf       - <baseresolution,Nspk,NSli>double - density compensation function
%           - num_spokes_chunk1 - int - (OPTIONAL) number spokes within the first scan session. 
%                                           If ignored, it is assumed to be equal to the number of spokes.
%
%       OUTPUTS:
%           
%           
%
%       DEPENDENCES:
%           - nufft_toolbox - ./external/alipour/utils/nufft_toolbox/ 
%                           citation:  Fessler, J. A., & Sutton, B. P. (2003). 
%                                       Nonuniform fast fourier transforms using min-max interpolation. 
%                                       IEEE Transactions on Signal Processing, 51(2), 560–574.
%           - MSMCNUFFT - ./external/@MCNUFFT
%                           by Jeff Fessler
%
%       AUTHOR:
%          Jaume Coll-Font <jcollfont@gmail.com>
%
%


function [rec_image_GRASP] = image_reconstruction_wrapper_GRASP(k_data, k_samples, dcf, coilprofile, spokes_per_vol, lambda, contrastSpokes)

    %% DEFINE
    [N, NCha, NSpk, NSli] = size(k_data);

    % create temporary folder
    tempFolder = tempname;
    mkdir(tempFolder);

    % eliminate the first spokes for contrast
    if ~exist('contrastSpokes')
        contrastSpokes = 1;
    end

    %% data selection
    % spokes_per_vol = 2*34;                         % number of lines (spokes) per volume
    valid_spoke_samples = 1:N;%floor(N/4+1):ceil(3*N/4);  % samples to use in each spoke    
    num_Volumes = floor( NSpk / spokes_per_vol); % number of volumes to be reconstructed
    valid_spokes = contrastSpokes:num_Volumes*spokes_per_vol;
    
    
    %% normalize coilprofile (redundant but does not hurt)
    maxCoilProfile = max(abs(coilprofile(:)));
    
    %% create cells for memory efficiency in parallel loop
    kdataCell = cell(1,NSli);
    for sl = 1:NSli
        kdataCell{sl} = 100*double(k_data(valid_spoke_samples,:,valid_spokes,sl));      % also make sure it is as double and rescale
    end
    clear k_data
    %%
    ksampCell = cell(1,NSli);
    for sl = 1:NSli
        ksampCell{sl} =  double(k_samples(valid_spoke_samples,valid_spokes,sl));
    end
    clear k_samples
    %%
    dcfCell = cell(1,NSli);
    for sl = 1:NSli
        dcfCell{sl} = double(dcf(valid_spoke_samples,valid_spokes,sl));
    end
    clear  dcf
    %%
    coilCell = cell(1,NSli);
    for sl = 1:NSli
        coilCell{sl} = double( coilprofile(:,:,:,sl) / maxCoilProfile );    % make sure it is double and normalize
    end
    clear coilprofile 

    %% Compute NUFFT matrix
    fprintf('\tGRASP for each slice:\n')
    rec_image_GRASP = cell(1,NSli);
    for sl = 1:NSli

        %% Solve inverse problem
        fprintf('\tComputing slice %d of %d\n',sl, NSli);
        [~, rec_image_GRASP{sl}, lambdaOut ] ...
                = function_grasp( ksampCell{sl}, coilCell{sl}, dcfCell{sl}, kdataCell{sl}, spokes_per_vol, lambda, true);

        % % free memory
        % ksampCell{sl} = [];
        % coilCell{sl} = [];
        % dcfCell{sl}= [];
        % kdataCell{sl} = [];
        % temporary save slice
        % tempGRASP_recon = rec_image_GRASP{sl};
        % tempSliceSave = sprintf( '%s/temp_GRASP_recon_slice%d.mat' ,tempFolder,sl);
        % save(  tempSliceSave, 'tempGRASP_recon' );

    end

    rec_image_GRASP = cell2mat( reshape(  rec_image_GRASP, [1,1,1,NSli] ) );
    
end
