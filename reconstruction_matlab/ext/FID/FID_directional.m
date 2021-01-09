%%
%
%      This function takes computes the auto FIDs from a stack of stars image and computes a weighted average to be sensitive to movement in the main coordinares
%
%
%       INPUT:
%           - k_data - - raw samples obtained from K-space
%           - coilprofile - -  precomputed coil profiles
%           - u - <3,D>double - each column is the desired sensitivity direction
%
function [FID_signal, W] = FID_directional(k_data, k_sample, dcf, coilprofile, u, targetArea)

    %% PARAMS
    [N, NCha, NSpk, NSli] = size(k_data);


    %% retrieve center of K-space in each spoke
    k_center_IX = round( N / 2 +1 ) + [-1:1];     % sample within spoke corresponding to center of K-space (assuming even # samples )
    T = 10;

    %% center FFT
    FID_signal = (fftshift(ifft(k_data(k_center_IX,:,:,:), NSli,1),1));

    %% introduce fourier transform
    fprintf('Preparing NUFFT\n')
    E = MCNUFFT( k_sample(k_center_IX,:,:), dcf(k_center_IX,:,:), coilprofile );
 
    %% compute weight vectors based on coil profiles
    fprintf('Computing weight matrix\n')
    [ W ] = compute_weightVectors( u, coilprofile, targetArea, E, k_center_IX, T, NSpk  );

    save('tmp/W2.mat','W','-v7.3');

    %% compute weighted sum of FID's
    fprintf('Computing FIDs\n')
    FID_signal = zeros( size(u,2), NSpk, NSli );
    for uu = 1:size(u,2)
        for sl = 1:NSli
            for tt = 1:NSpk
                current_spokes = tt:min(NSpk,tt+T-1);
                FID_signal(uu,tt,sl) = abs( sum(sum(sum(   squeeze(W(uu,:,sl,:,1:numel(current_spokes))) .* squeeze(k_data(:,:, current_spokes, sl)) ) )));
            end
        end
    end

end

%%
%
%   This sub-function computes the weight vectors used to maximize the FID sensitivity to the desired directions u.
%
%   INPUT:
%       - u                 - <3,D>double - each column is the desired sensitivity direction
%       - coilprofile       -           - coil profiles.
%       - target area       - <>        - area of the image to be targeted. (1/0)  mask
%       - E                 - MCNUFFT class - class computing the NUFFT for the multicoil image.
%       - valid_K           - <1,M> int - indices of the valid samples in each spoke.
%       - T                 - int       - number of spokes to use in each window.
%
function [ W ] = compute_weightVectors( u, coilprofile ,targetArea, E, valid_K, T, NSpk )

    numRep = 100;
    lambda = 1;

    alpha_int = 1.1;    % increment in intensity
    vx_step = 1;       % size of voxel change 

    % define 
    [Nx,Ny,NCha,NSli] = size(coilprofile);
    M = numel(valid_K);

    % vectorize and normalize coil
    coil_v =  permute(coilprofile,[1,2,4,3]);
    coil_v = coil_v ./ repmat( sqrt( sum(sum(sum( abs(coil_v).^2 ,1),2),3) ) ,[Nx,Ny,NSli,1]);


    W = zeros( NCha*T*M*NSli ,size(u,2));
    for dd = 1:size(u,2)
        
        fprintf('\tDirection %d of %d\n', dd, size(u,2))
        %% create target images whose change corresponds to the desired direction
        % initialize target images with random distribution of points
        R = zeros(NCha*T*M*NSli, NCha*T*M*NSli);
        for  rr = 1:numRep

            p = round(diag([Nx,Ny,NSli])*rand(3,floor(Nx*Ny*NSli/2)));
            
            S = cell(1,T);
            vR = zeros(0,NSpk-T);
            vStat = zeros(0,NSpk-T);
            for tt = 1:T
                
                timeSamp = [1:(NSpk-T)]  + tt;

                % prepare target images
                It = zeros(Nx,Ny,NSli);
                for ii = 1:size(p,2);               % for every voxel
                    dp = p(:,ii) + vx_step*(tt-1)*u(:,dd);
                    dp = min( [Nx;Ny;NSli]  , dp );
                    dp = max( [1;1;1]       , dp );
                    It(dp(1),dp(2),dp(3)) = 1;
                end
                It = It.*targetArea;
                
                %% generate "correlation" matrix to maximize the directivity
                tempS = permute(    E * It,   [1,4,3,2]); % [N, NCha, NSli, NSpk]
                S{tt} = reshape(  tempS(:,:,:,timeSamp), [M*NSli*NCha,NSpk-T]);  % [M *NCha*NSli, NSpk]
                vR = [ vR; S{tt}];  % [M *NCha *NSli  * T, NSpk]

                %% generate "correlation" matrix to minimize static change in amplitude
                vStat = [ vStat; alpha_int^(tt)*S{1}]; 
            end

            %% add regularization to avoid being sensitive to the image in itself
            R = R + vR*vR'/numRep  - lambda/numRep*vStat*vStat';  % [M *NCha *NSli  * T, M *NCha *NSli  * T]

        end

        save(sprintf('tmp/R%d.mat', dd),'R','-v7.3');

        %% solve maximization problem (take the biggest eigenvector)
        [v,s] = eigs(R);
        W(:,dd) = v(:,1);

    end
    % reshape for correct output size
    W = permute(  reshape( W, [ M, NSli, NCha, T, size(u,2)] ), [5,1,2,3,4]);    % [Dirs, M, NSli, NCha, T ]

end


%%
%
%
%
%
%
% function [] = compute_NUFFT()
% disp 'gato';
% end