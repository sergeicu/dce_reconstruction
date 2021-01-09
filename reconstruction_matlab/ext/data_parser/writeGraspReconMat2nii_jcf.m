%% writeGraspRecon2mat2nii(patientName0)
% this function reads the saved grasp resonstructed slices writes it to .nii
% files
%
%
%       params:
%           - timePoints 
%           - dirNameOut
%           - saveName
%           - keepSlices - DO NOT ADD IT!

%
%%
function []=writeGraspReconMat2nii_jcf( data, patientName0,dirNameOut, saveName,  sliceOverSampling, referenceNii, indivVolumes)

    % dirNameOut='/fileserver/abd/GraspRecons/reconResultsGRASP/method2/';

    %% read headers from scanner reconstructed data
    a1 = load_untouch_nii(referenceNii);
    dim=a1.hdr.dime.dim;
    res=a1.hdr.dime.pixdim;



    %% reconstruct vol4D from saved slices
    tdim = numel(data);
    [xdim,ydim,zdim]=size( data{1} );
    startxy = ((xdim-dim(2))/2)+1;
    endxy = ((xdim-dim(2))/2)+dim(2);

    temp = cell2mat(reshape(data,[1,1,1,tdim]));
    if xdim ~= 448
        temp = imresize( temp, [448, 448] );
    end
    
    [xdim,ydim,zdim,tdim]=size( temp );
    startxy = ((xdim-dim(2))/2)+1;
    endxy = ((xdim-dim(2))/2)+dim(2);

    vol4D = zeros(dim(2),dim(3), zdim, tdim);
    vol4D = temp(startxy:endxy,startxy:endxy,:,:);
   
    clear temp;

    %% fft shift in z dimention (for each volume)
    for tNdx=1:tdim
        vol4D(:,:,:,tNdx)=abs(fftshift(vol4D(:,:,:,tNdx),3));
    end

    
    %% adjust z-dim resolution based on scanner recon resolution
    totalslices=round(dim(4)*(1+sliceOverSampling));    
    dimDiff=round(dim(4)*(1+sliceOverSampling))-zdim; % sila: Consider slice oversampling
    if mod(dimDiff,2)==0
        ZpadBoth=dimDiff/2;
        vol4D2 = abs(ifft(ifftshift(padarray(fftshift(fft(vol4D,[],3),3),[0 0 ZpadBoth 0]),3),[],3));
    elseif mod(dimDiff,2)==1
        Zpad1=floor(dimDiff/2);Zpad2=dimDiff-Zpad1;
        padded1=padarray(fftshift(fft(vol4D,[],3),3),[0 0 Zpad1 0],'pre');
        padded=padarray(padded1,[0 0 Zpad2 0],'post');
        vol4D2 = abs(ifft(ifftshift(padded,3),[],3));
    end

    vol4D3=vol4D2;  % added by Sila
    vol4D2=vol4D3(:,:,floor((totalslices-dim(4))/2)+1:dim(4)+floor((totalslices-dim(4))/2),:);     % added by Sila    
 
    %%%% rescale the image for uint16
    maxPixels=max(vol4D2(:));
    minPixels=min(vol4D2(:));
    vol4D = round( (vol4D2 -minPixels) ./ (maxPixels-minPixels) * 2^16 );

    % minPixels=min(vol4D2(:));maxPixels=max(vol4D2(:));
    % m = 2^16/(max(maxPixels) - minPixels);
    % vol4D = imlincomb(m, vol4D2, -(m * minPixels), 'uint16');
    
    %% save .nii file for each volume    
%     if ~exist('indivVolumes')
        if indivVolumes
            for tNdx=1:tdim
                a1.img=vol4D(:,:,:,tNdx);
                save_untouch_nii(a1,strcat(dirNameOut,'/',saveName, '_',num2str(tNdx-1),'.nii'));
            end
        end
%     end
    a1.img=vol4D;
    dim(5)=size(vol4D,4);
    dim(1)=4;
    a1.hdr.dime.dim=dim;
    res(5)=1;
    a1.hdr.dime.pixdim=res;
    save_untouch_nii(a1,strcat(dirNameOut,'/',saveName,'.nii'));
    
end