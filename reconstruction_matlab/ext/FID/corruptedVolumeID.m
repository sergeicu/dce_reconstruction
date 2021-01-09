function [outlierMaskVol, outlierMaskFID] = corruptedVolumeID(  corrFID, spk_per_vol, T, startIntakeVols )

    corrFID = max(1-corrFID,[],1);
    
    NSpk = numel(corrFID);
    
    corrFID_thr =  max( 0.01,min( 0.1, std(corrFID(:)))  );    
    
    movementIX = find(corrFID > corrFID_thr);
    outlierMaskFID = zeros(1,NSpk);
    outlierMaskFID(movementIX) = 1;
    
    % eliminate the intake period as possible outlier
    outlierMaskFID(1:spk_per_vol*startIntakeVols) = 0;

    % filter data with opening filter
	SE = strel('disk',5);
    outlierMaskFID = imclose( outlierMaskFID, SE );

    
    outlierMaskVol = zeros(1,T);
    for rr =1:T
        aa = [spk_per_vol*(rr-1)+1:spk_per_vol*rr];
        if sum( outlierMaskFID(aa) ) > 0
            outlierMaskVol(rr) = 1;
        end
    end


end