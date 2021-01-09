%% 
%
%       This function applies the reconstruction of the MR image from 
%       k-space sampling. The method used for reconstruction is called GRASP:
%       
%           L. Feng, R. Grimm, K. T. obia. Block, H. Chandarana, 
%           S. Kim, J. Xu, L. Axel, D. K. Sodickson, and R. Otazo, 
%           “Golden-angle radial sparse parallel MRI: combination of compressed sensing, parallel imaging, 
%           and golden-angle radial sampling for fast and flexible dynamic volumetric MRI,” 
%           Magn. Reson. Med., vol. 72, no. 3, pp. 707–717, 2014.
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


function [rec_image] = image_reconstruction_GRASP( k_data, k_samples_1, dcf_1, coilprofile, num_spokes_chunk1 )

    %% DEFINE
    [N, NCha, NSpk, NSli] = size(k_data);

    %% data selection
    spokes_per_vol = 34;                         % number of lines (spokes) per volume
    valid_spoke_samples = (N/4+1):(3*N/4);  % samples to use in each spoke    
    num_Volumes = floor( NSpk / spokes_per_vol); % number of volumes to be reconstructed


    k_data_perm = permute(k_data,[1,3,2,4]); % permuting the data to  samples / spokes / chanels / slices

    %% for each volume Solve inverse problem  without regularization (initialize problem)
    fprintf('Running Standard Inverse Problem on Volumes:\n')
    rec_image = cell(1,num_Volumes);
    for vv = 1:num_Volumes

        %% select data to use
        current_spokes  = (vv-1)*spokes_per_vol+1:min(vv*spokes_per_vol ,NSpk);                 % spokes to use for current volume
  
        %% Compute NUFFT matrix
        fprintf('\tNUFFT for volume: %d out of %d\n', vv, num_Volumes)
        E = MCNUFFT( double(k_samples( valid_spoke_samples, current_spokes, : )),...
                     double(dcf( valid_spoke_samples, current_spokes, : )), ...
                     double(coilprofile(valid_spoke_samples,valid_spoke_samples,:,:)) );

        %% Solve inverse problem
        fprintf('\tComputing image\n')
        rec_image{vv} = E'*k_data_perm( valid_spoke_samples, current_spokes, :, :);

    end

    %% run GRASP on all volumes

end


%%
%
%   This subfunction implements a Total Variation optimization using on Conjugate Gradient.
%
%
function [] = conjugateGradient_TVoptimization()



end