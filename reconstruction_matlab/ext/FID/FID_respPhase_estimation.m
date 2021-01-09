%%
%
%       This function estimates the respiratory phase from the FID signal.
%
%
function [Res_Signal] = FID_respPhase_estimation( data )

    t = (size(data,1)/2+1);
    center_profile = zeros([size(data,2) size(data,3) size(data,4)]);
    center_profile(:,:,:) = data(t,:,:,:);
    center_profile=permute(center_profile,[3 2 1]);
    kc = center_profile;

    nc=30;
    ZIP=abs(fftshift(ifft(kc,400,1),1));
    [nz,ntviews,nc]=size(ZIP);
    for ii=1:nc
        for jj=1:ntviews
            maxprof=max(ZIP(:,jj,ii));
            minprof=min(ZIP(:,jj,ii));
            ZIP(:,jj,ii)=(ZIP(:,jj,ii)-minprof)./(maxprof-minprof);
        end
    end


     ZIP1=ZIP(:,1:end,:);
    [nz,ntviews,nc]=size(ZIP1);

    % Perform PCA on each coil element
    kk=1;clear PCs
    for ii=1:nc
        tmp=permute(ZIP1(:,:,ii),[1,3,2]);
        tmp=abs(reshape(tmp,[size(tmp,1)*size(tmp,2),ntviews])');

        covariance=cov(tmp);
        [tmp2, V] = eig(covariance);
        V = diag(V);
        [junk, rindices] = sort(-1*V);
        V = V(rindices);
        tmp2 = tmp2(:,rindices);
        PC = (tmp2' * tmp')';

        % Take the first two principal components from each coil element.
        for jj=1:2
            tmp3=smooth(PC(:,jj),6,'lowess');
            tmp3=tmp3-min(tmp3(:));
            tmp3=tmp3./max(tmp3(:));
            PCs(:,kk)=tmp3;kk=kk+1;

        end
    end    



    %% correlation coefficient

    d1=PCs;
    for i = 2 : length(PCs)
             cc=corrcoef(d1(i-1,:)',d1(i,:)');
            corrt2(i) = cc(1,2);
    end
    corrt2(1)=corrt2(2);
    
    %%
    thresh = 0.95;
    [Res_Signal, cluster] = CoilClustering(PCs, thresh);
  

end