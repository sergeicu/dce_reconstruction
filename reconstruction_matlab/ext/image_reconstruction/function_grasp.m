% GRASP (Golden-angle RAdial Sparse Parallel MRI)
% Combination of compressed sensing, parallel imaging and golden-angle
% radial sampling for fast dynamic MRI.
% This demo will reconstruct one slice of a contrast-enhanced liver scan.
% Radial k-space data are continously acquired using the golden-angle
% scheme and retrospectively sorted into a time-series of images using
% a number of consecutive spokes to form each frame. In this example, 21
% consecutive spokes are used for each frame, which provides 28 temporal
% frames. The reconstructed image matrix for each frame is 384x384. The
% undersampling factor is 384/21=18.2. TheNUFFT toolbox from Jeff Fessler
% is employed to reconstruct radial data
%
% Li Feng, Ricardo Otazo, NYU, 2012
function [recon_nufft,recon_cs,lambdaOut]=function_grasp(k, coilprofile, dcf, data,nspokes,lambda,regEnabled)
    % clear all;close all;[k, coilprofile, dcf, data]
    % addpath('nufft_toolbox/');
    % display('hello!')
    %load siladatax;load siladcf;load silacoilprofile;load silak
    % define number of spokes to be used per frame (Fibonacci number)
    %nspokes=31;
    % load radial data
    % load liver_data.mat
    % b1=b1/max(abs(b1(:)));
    % data dimensions
    %kdata=double(datax(1:2:end,:,:));
    % kdata=double(data(:,:,:));
    % clear data;
    kdatax=permute(double(data),[1 3 2]);       % Data dimensions:  [ num points per line  |  number of lines  | number of coil channels ]
    [nx,ntviews,nc]=size(kdatax);
    %w=ones(size(kdatax(:,:,1)));
    w=dcf;
    for ch=1:nc,kdatax(:,:,ch)=kdatax(:,:,ch).*sqrt(w);end      
    % number of frames
    nt=floor(ntviews/nspokes);
    % crop the data according to the number of spokes per frame
    kdata1=kdatax(:,1:nt*nspokes,:);
    clear kdatax;
    k=k(:,1:nt*nspokes);
    w=w(:,1:nt*nspokes);

    %b1=coilprofile(113:112+224,113:112+224,:);
    % b1=coilprofile;
    % b1=ones(224,224,34);
    %b1=abs(b1);
    coilprofile=coilprofile/max(abs(coilprofile(:)));
    %b1=coilprofile;
    % sort the data into a time-series
    kdatau=zeros(size(kdata1,1),nspokes,size(kdata1,3));
    ku=zeros(size(k,1),nspokes);
    wu=zeros(size(w,1),nspokes);
    for ii=1:nt
        kdatau(:,:,:,ii)=kdata1(:,(ii-1)*nspokes+1:ii*nspokes,:);
        ku(:,:,ii)=k(:,(ii-1)*nspokes+1:ii*nspokes);
        wu(:,:,ii)=w(:,(ii-1)*nspokes+1:ii*nspokes);
    end

    clearvars k w kdata1;
    % multicoil NUFFT operator
    param.E=MCNUFFT(ku,wu,coilprofile);
    % undersampled data
    param.y=double(kdatau);
    %clear kdata kdatau k ku wu w
    % nufft recon
    recon_nufft=param.E'*double(param.y);
%     recon_nufft = zeros(size(recon_nufft));
    % param.y=double(kdatau(:,:,1));
    % param.E=MCNUFFT(ku,wu,ones(448,448));
    %
    % recon_nufft2=param.E'*double(param.y);
    % % parameters for reconstruction
    display('hello!')
    if regEnabled
    %     display('hello!')
        param.W = TV_Temp();
        
        if isempty(lambda)
            param.lambda=0.0125*max(abs(recon_nufft(:)));
        else
            param.lambda=lambda*max(abs(recon_nufft(:)));
        end
        lambdaOut=param.lambda;
        
        param.nite = 8;
        param.display=1;
        fprintf('\n GRASP reconstruction \n')
        
        recon_cs=recon_nufft;
        for n=1:3
            recon_cs = CSL1NlCg(recon_cs,param);
        end    
        
    else
        recon_cs=recon_nufft;
        lambdaOut=0;
    %     display('hello!')
    end
    % [user, ~] = memory;
    % display(user.MemUsedMATLAB);
    % recon_nufft=flipdim(recon_nufft,1);
    % recon_cs=flipdim(recon_cs,1);
    %
    % % display 4 frames
    % recon_nufft2=recon_nufft(:,:,1);recon_nufft2=cat(2,recon_nufft2,recon_nufft(:,:,7));recon_nufft2=cat(2,recon_nufft2,recon_nufft(:,:,13));recon_nufft2=cat(2,recon_nufft2,recon_nufft(:,:,23));
    % recon_cs2=recon_cs(:,:,1);recon_cs2=cat(2,recon_cs2,recon_cs(:,:,7));recon_sizecs2=cat(2,recon_cs2,recon_cs(:,:,13));recon_cs2=cat(2,recon_cs2,recon_cs(:,:,23));
    % figure;
    % subplot(2,1,1),imshow(abs(recon_nufft2),[]);title('Zero-filled FFT')
    % subplot(2,1,2),imshow(abs(recon_cs2),[]);title('GRASP')
end