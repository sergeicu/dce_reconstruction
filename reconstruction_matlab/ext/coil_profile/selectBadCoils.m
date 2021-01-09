%% HELP:
%
%       This function selects which coils are corrupted with artifacts.
%       It does this task by evaluating "how diagonal" is the local noise covariance matrix of each coil reconstructed image.
%       
%       The algorithm first computes an estimate of the covariance of a patch on the reconstructed images of size covDim. 
%       This covariance is estimated as the moving average covariance along a wider patch of the image edges of size W. (assming these are background pixels)
%       Then, the algorithm evaluates the matrix 2-norm of the estimated covariance matrix minus a diagonal. 
%       (the larger the difference the less diagonal-like the matrix is)
%
%       INPUT:
%           - coil_rec_image - <Nx,Ny,Nz,NCha>complex - the reconstruted coil images.
%           - covDim - <1,2>int - size of the local patch to compute covariance matrix (e.g. [10,10]).
%           - W - <1,2>int - background pixels window taken at the edges (e.g. [100,100]).
%           - err_th - double - error threshold to use.
%           
%       OUTPUT:
%           - goodCoils - <1,NGoodCoils>int - coils with "non-diagonal-like" metric below threshold.
%           - covErr - <1,NCha>double - "non-diagonal-like" metric .
%
%       AUTHOR:
%           Jaume Coll-Font <jaume.coll-font@childrens.harvard.edu>
%
%%
function [goodCoils, covErr] = selectBadCoils( coil_rec_image, covDim, W, err_th )

    
    [Nx, Ny, NSli, NCha] = size(coil_rec_image);
    
    coilCell = cell(1,NCha);
    for ch = 1:NCha
        coilCell{ch} = coil_rec_image(:,:,:,ch);
    end
    clear 'coil_rec_image';

    covErr = zeros(NSli,NCha);
    goodCoils = ones(NSli, NCha);
    parfor ch =1 :NCha

        
        for sl = 1:NSli
            movCov = zeros(covDim);
            for edg = 1:4       % for all 4 edges
                switch edg
                    case 1
                        r_conv = [floor(covDim(1)/2)+1:W(1)-floor(covDim(1)/2)];
                        c_conv = [floor(covDim(2)/2)+1:W(2)-floor(covDim(2)/2)];
                    case 2
                        r_conv = [floor(covDim(1)/2)+1:W(1)-floor(covDim(1)/2)];
                        c_conv = Ny - W(2) + [floor(covDim(2)/2)+1:W(2)-floor(covDim(2)/2)];
                    case 3
                        r_conv = Nx - W(1) + [floor(covDim(1)/2)+1:W(1)-floor(covDim(1)/2)];
                        c_conv = [floor(covDim(2)/2)+1:W(2)-floor(covDim(2)/2)];
                    case 4
                        r_conv = Nx - W(1) + [floor(covDim(1)/2)+1:W(1)-floor(covDim(1)/2)];
                        c_conv = Ny - W(2) + [floor(covDim(2)/2)+1:W(2)-floor(covDim(2)/2)];
                end
                for r = r_conv
                    for c = c_conv
                        I = coilCell{ch}(r + [-floor(covDim(1)/2):floor(covDim(1)/2)-1],c + [-floor(covDim(2)/2):floor(covDim(2)/2)-1], sl );
                        movCov = movCov + corr(I);
                    end
                end
            end

            % normalize covariance
            movCov = movCov / mean(diag(movCov));

            covErr(sl,ch) = norm( movCov - eye(covDim) ,2);
        end        
    end

    goodCoils = mean(abs(covErr),1) < err_th;

end