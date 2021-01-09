%%
%   FID_extraction
%       
%       This function takes data sampled as a stack of spokes and computes an FID movement tracking signal.
%
%          The main assumption is that the mid-samples in each spoke correspond to the center in K-space. 
%          The change of value of these centers indicate the variation in position of the subject.
%
%       INPUTS:
%           - k_data - <N,NCha,NSpk,NSli>double - data sampled in K-space.
%                                               Dimensions: - N = number of samples per spoke.
%                                                           - NCha = number of chanels.
%                                                           - NSpk = number of spokes.
%                                                           - NSli = number of slices.
%
%
%       OUTPUTS:
%           - FID_signal   - <3,T>double - 
%           - FID_corr      - <>double - 
%       
%
function [FID_signal, FID_corr, FID_meanNormChange, FID_maxNormChange, FID_medNormChange] = FID_extraction( k_data )

    %% PARAMS
    [N, NCha, NSpk, NSli] = size(k_data);
    
    %% retrieve center of K-space in each spoke
    k_center_IX = round( N / 2 +1 );     % sample within spoke corresponding to center of K-space (assuming even # samples )
    k_center = permute( squeeze( k_data( k_center_IX, :,:,: ) ) ,[3,2,1]); % [Nsli,Nspk,NCha]


    %% Interpolate data to extended number of slices and normalize between 0 and 1
    % FID_signal = abs(k_center);
    FID_signal=abs(fftshift(ifft(k_center,NSli,1),1));      % interpolate data to new location
    maxprof = squeeze( max(FID_signal,[],1) ); 
    minprof = squeeze( min(FID_signal,[],1) );
    
    for ii=1:NCha
        for jj=1:NSpk
            FID_signal(:,jj,ii)=(FID_signal(:,jj,ii)-minprof(jj,ii))./(maxprof(jj,ii)-minprof(jj,ii));
        end
    end


    % filter with a moving mean average
    meanFID = movmean(abs(squeeze(FID_signal)),34,2);

    % compute correlation FID
    step = 17;
    FID_corr = zeros( size(FID_signal,1),size(FID_signal,2));
    for sl = 1:size(FID_signal,1)
        for ii = 1:size(FID_signal,2)
            FID_corr(sl,ii) = -1;           % -1 compensates for the correlation with self, which would bias results
            for w = max(1,ii-step):min(size(FID_signal,2),ii+step)            % displace along a window
                cc = corrcoef(meanFID(sl,ii,:), meanFID(sl,w,:));
                FID_corr(sl,ii) = FID_corr(sl,ii) + cc(1,2);
            end
            FID_corr(sl,ii) = FID_corr(sl,ii) / (min(size(FID_signal,2),ii+step)-max(0,ii-step));
        end
    end

    % compute time change change
    FIDchange = diff( meanFID ,1,2);
    FIDchange = FIDchange ./ meanFID(:,1:end-1,:);
    % compute normalized mean change
    FID_meanNormChange =  mean( FIDchange ,3);
    % compute max normalized change
    FID_maxNormChange =  max( FIDchange ,[],3);
    % compute median real normalized change
    FID_medNormChange = median( abs( FIDchange ) ,3);

end