%%
%
%       This function loads the K-space data from a .dat file and sets it up in the appropriate form.
%       The code assumes that the loaded data has been captured as a stack of stars (stokes) with 
%       multiple coils and possibly over a series (usually 2) of recordings.
%
%
%       This function is copied from Sila's GRASP script as is.
%
%       INPUT:
%           - file - <1,F>cell - each cell contains the path <string> to the data of each scan session.
%           - timePoints - <1,F>int - time points per scan session (num of volumes) assuming a fixed (yet unknown) frame rate.
%
%       OUTPUT:
%           - k22n - <N,NCha,NSpk,NSli>double - data sampled in K-space.
%                                               Dimensions: - N = number of samples per spoke.
%                                                           - NCha = number of chanels.
%                                                           - NSpk = number of spokes.
%                                                           - NSli = number of slices.
%           - k3n - <N,NCha,NSpk,NSli>double - data sampled in K-space Fourier transformed along slice dimension.
%                                               Dimensions: - N = number of samples per spoke.
%                                                           - NCha = number of chanels.
%                                                           - NSpk = number of spokes.
%                                                           - NSli = number of slices.
%           - k_samples - <baseresolution,1>double -    k space coordinates ascomplex values (kx + i ky)
%           - dcf       - <baseresolution,NSpk,NSli>double - density compensation function
%           - totalSpokesInSeq - cell - number of spokes in each sequence
%
%
function [  k22n, k3n, k_samples, dcf, totalSpokesInSeq, SliceOversampling, coilList2] = load_kspace_data( file, timePoints, spokesPerVol  )


    %% JCF: from readParamsGrasp

    %% read and cat data
    % read each raw file and concatenate all chuncks in a raw file (part 1 or 2) 
    for i=1:size(file,1)
        k2=mapVBVD(char(file(i,:)));

        if iscell(k2)
            k2 = k2{2};
        end
        if size(k2,2)>1   %if it is binney7
            k2=k2{2};  
        end
        
        SliceOversampling = k2.hdr.MeasYaps.sKSpace.dSliceOversamplingForDialog;
        
        
        k22=k2.image{''};
        scannerModel=k2.hdr.Dicom.ManufacturersModelName;
        dicomfolderName=k2.hdr.Config.ProtocolName;
        
        if strcmp(scannerModel,'TrioTim') || strcmp(scannerModel,'Avanto')
            for j=1:size(k2.hdr.MeasYaps.asCoilSelectMeas{1,1}.asList,2)
                coilList{i,j}=k2.hdr.MeasYaps.asCoilSelectMeas{1,1}.asList{1,j}.sCoilElementID.tElement(2:end-1);
            end
        elseif strcmp(scannerModel,'Skyra') || strcmp(scannerModel,'Prisma')
            for j=1:size(k2.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.asList,2)
                coilList{i,j}=k2.hdr.MeasYaps.sCoilSelectMeas.aRxCoilSelectData{1,1}.asList{1,j}.sCoilElementID.tElement(2:end-1);
            end
        end
        
        
        numOfCoils=size(k22,2);
        numOfSpokes=size(k22,3);
        numofSlices=size(k22,4);
        numOfChunks=size(k22,5);
        
        % chunk concatenation
        k22n=zeros(size(k22,1),size(k22,2),numOfSpokes*numOfChunks,numofSlices);
        for j=1:numOfChunks
            k22n(:,:,((j-1)*numOfSpokes+1):j*numOfSpokes,:)=k22(:,:,:,:,j);
        end
        
        
        k22nall{i}=k22n;
        numOfSpokesA(i)=numOfSpokes;
        numOfChunksA(i)=numOfChunks;
        numOfCoilsA(i)=numOfCoils;
        totalSpokesInSeq{i}=repmat(numOfSpokes,[1,numOfChunks]);
        clear k2 k22 k22n;
    end

    %%
    if exist('spokesPerVol')
        if numel(spokesPerVol) ==0
            spokesPerVol=floor(size(k22nall{1},3)/timePoints(1));
        end
    end

    if length(unique(numOfCoilsA))>1
        display('COIL MISMATCH!!!');
    %     k22nall2=coilMismatchCorrection(k22nall);
        [nOfCoils,ndx]=min(numOfCoilsA);
        selectedCoils=coilList(ndx,1:nOfCoils);
        coilList2(ndx,:)= selectedCoils;
        extraCoilSeqs=setdiff(1:length(numOfCoilsA),ndx);
        for ci=extraCoilSeqs
            commonNdc=ismember(coilList(extraCoilSeqs,:),selectedCoils);
            k22nall{ci}=k22nall{ci}(:,commonNdc,:,:);
    %         coilListTemp=coilList(ci,commonNdc);
            coilList2(ci,:)= coilList(ci,commonNdc);
        end
        
    %     if ~isequal()
    else
        coilList2=coilList;
    end


    if size(k22nall,2)>1
        %     concatenateThem;
        k22n=k22nall{1};
        %%% check for coil order match 
        firstSeqOrder=cellstr(coilList2(1,:));
        for i=2:size(k22nall,2)
            seqCoilList=coilList2(i,:);
    %         [~,ndxx]=ismember(firstSeqOrder,cellstr(seqCoilList(~cellfun('isempty',seqCoilList))));
            [~,ndxx]=ismember(firstSeqOrder,cellstr(seqCoilList));
            if ~isequal(ndxx,sort(ndxx))
            k22nall{i}=k22nall{i}(:,ndxx,:,:);
            end        
            spokesselect=floor(size(k22nall{i},3)/spokesPerVol)*spokesPerVol;
            if ~isequal(size(k22nall{1},1),size(k22nall{i},1))
                padBoth=(size(k22nall{1},1)-size(k22nall{i},1))/2;
                k22nall{i} = padarray(k22nall{i},[padBoth 0 0 0]);
            end
            k22n=cat(3,k22n,k22nall{i}(:,:,1:spokesselect,:));
            totalSpokesInSeq{i}(numOfChunksA(i))=spokesselect-(numOfChunksA(i)-1)*numOfSpokesA(i);
        end
    else
        k22n=k22nall{1};
    end
    clear k22nall;
    %%
    % [dim1: number of points on each radial line(spoke)] [dim2: number of channels]
    % [dim3: number of spokes] [dim4: number of slices(volumes)] [dim5: number of chunks]
    % number slices


    %%  save slices to file
    k3n=fft(double(k22n),[],4);


    %% positions in K-space
    [N, NCha, NSpk, NSli] = size(k3n);


    k_samples=zeros(N,NSpk, NSli);
    dcf=zeros(N,NSpk, NSli);
    for sl = 1:NSli
        sp1=[0 cell2mat(totalSpokesInSeq)];
        sp2=cumsum(sp1);
        for i=1:length(sp1)-1
            [ k_samples(:,sp2(i)+1:sp2(i+1),sl), dcf(:,sp2(i)+1:sp2(i+1),sl)] =  calc_k_dcf(100*k3n(:,:,sp2(i)+1:sp2(i+1),sl),1);
        end
    end

end