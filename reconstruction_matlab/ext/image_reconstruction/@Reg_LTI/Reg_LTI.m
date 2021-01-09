function  res = Reg_LTI( generatedLTIDataPath, genTag, maskPath, sliceImg, sliceK )

    % define the LTI loss function
    
    % if no tag, no worries
    if ~exist('genTag')
        genTag = '.';
    end
    
    % load data generated using LTI
%     [genFilesList] = loadGenFiles( generatedLTIDataPath, genTag );
%     [genLTIData] = loadGeneratedData(genFilesList);
%     genLTIData = load_untouch_nii(generatedLTIDataPath);
      
    load(generatedLTIDataPath);
    genLTIData = ltifit;
    
    % adapt mask size to image size  (HARDCODED AND ASSUMING THAT IMAGES
    % ARE CROPPED DURING SAVING PROCEDURE) .... sorry future me
    mask = load_untouch_nii(maskPath{1});
    bigMask = mask.img(:,:,sliceImg);
    % for mm = 2:numel(maskPath)
    %     mask = load_untouch_nii(maskPath{mm});
    %     bigMask = bigMask + mask.img(:,:,sliceImg);
    % end
    bigMask = ones(size(mask.img(:,:,sliceImg)));
    
    % set parameters
    [Nx, Ny, Nz, T] = size(genLTIData);
    res.resol = [Nx, Ny, Nz, T];
    res.LTIData = double(reshape(squeeze(genLTIData(:,:,sliceK,:)),[Nx*Ny*numel(sliceImg),T]));       % manua normalization, should be based on the NUFFT
    res.maskIX = find(bigMask(:));
    res.adjoint = 0;
    res = class(res,'Reg_LTI');

end


function [genFilesList] = loadGenFiles( generatedLTIDataPath, genTag )

    allFiles = dir(generatedLTIDataPath);
    
    genFilesList = {};
    for ii = 1:numel(allFiles)
        if ( strfind(allFiles(ii).name, '.nrrd') > 0 ) & ( strfind(allFiles(ii).name, genTag) > 0 )
            genFilesList = [genFilesList , { [allFiles(ii).folder, '/' , allFiles(ii).name] }];
        end
    end
end


function [genLTIData] = loadGeneratedData(genFilesList)

    Nfiles = numel(genFilesList);
    
    for ii = 1:Nfiles
        tmpData = nrrdread( genFilesList{ii} );
        tmpData = permute(tmpData,[2,3,1]);         % X, Y, Z dimensions (nrrdread loads Z first)
        tmpData = tmpData(:,end:-1:1,:);         % Y dimension is saved in the opposite directionenvinGuda1
        
        if ii == 1
            [X,Y,Z] = size(tmpData);
            genLTIData = zeros(X,Y,Z,Nfiles);
        end
        genLTIData(:,:,:,ii) = tmpData;
    end

end
