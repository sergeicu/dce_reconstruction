
function [recon_nufft,recon_cs,lambdaOut]=function_L1LTI( inputDataStruct, nspokes, lambda, regEnabled, regularizationClass)


    %% retrieve basic params
    [nx,nc,ntviews]=size(inputDataStruct.data);       % dimensions of data
    nt=floor(ntviews/nspokes);          % number of frames
    
    
    %% RESHAPE DATA INTO APPROPROATE WAY (permute and crop over line dimension)
    kdatax = permute(double(inputDataStruct.data(:,:,1:nt*nspokes)) ,[1 3 2]); % Data dimensions:  [ num points per line  |  number of lines (cropped to have integer num frames)  | number of coil channels ]
    w = inputDataStruct.dcf(:,1:nt*nspokes,:);
    k = inputDataStruct.k(:,1:nt*nspokes);      
    
    % pre-weight data with density compensation function
    for ch=1:nc
        kdatax(:,:,ch) = kdatax(:,:,ch).*sqrt(w);
    end    
    
    % normalize coil profiles
    coilprofile = inputDataStruct.coilprofile/max(abs(inputDataStruct.coilprofile(:)));


    % sort the data into a time-series
    kdatau = zeros( nx, nspokes, nc);
    ku = zeros( nx,nspokes);
    wu = zeros( nx,nspokes);
    for ii=1:nt
        kdatau(:,:,:,ii)=kdatax(:,(ii-1)*nspokes+1:ii*nspokes,:);
        ku(:,:,ii)=k(:,(ii-1)*nspokes+1:ii*nspokes);
        wu(:,:,ii)=w(:,(ii-1)*nspokes+1:ii*nspokes);
    end

    clearvars k w kdatax;
    
    %% CREATE NUFFT OPERATOR AND SAVE EVERYTHING INTO STRUCT MODE
    % multicoil NUFFT operator
    param.E = MCNUFFT(ku,wu,coilprofile);
    
    % undersampled data
    param.y = double(kdatau);
    
    clear kdatau ku wu
    
   
    %% INITIALIZE WITH NUFFT RECONSTRUCTION
    recon_nufft = param.E'*double(param.y);
 
    %% RUN!!!
    display('hello!')
    if regEnabled
        
%         param.W = Reg_LTI( regularizationStruct.generatedLTIDataPath, regularizationStruct.genTag, regularizationStruct.maskPath, regularizationStruct.slice);
        param.W = TV_Temp();
        param.W2 = regularizationClass;
        
        if isempty(lambda)
            param.lambda=0.0125*max(abs(recon_nufft(:)));
            param.lambda2=0.0125*max(abs(recon_nufft(:)));
        else
            param.lambda=lambda(1)*max(abs(recon_nufft(:)));
            param.lambda2=lambda(2)*max(abs(recon_nufft(:)));
        end
        lambdaOut=[param.lambda,param.lambda2];
        
        param.nite = 8;
        param.display=1;
        fprintf('\n LTI based reconstruction \n')
        
        recon_cs=recon_nufft;
        for n=1:3
            recon_cs = CSLTINlCg(recon_cs,param);
        end    
        
    else
        recon_cs=recon_nufft;
        lambdaOut=0;
    end
end