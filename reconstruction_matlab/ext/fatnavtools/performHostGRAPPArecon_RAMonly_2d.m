function [grappaRecon, mOutGRAPPA] = performHostGRAPPArecon_RAMonly_2d(twix_obj, iRepToRecon, iEcoToRecon, nSliceNeighbours)
%
% function [grappaRecon_1DFFT, mOutGRAPPA] = performHostGRAPPArecon_RAMonly(twix_obj)
%
% Applies GRAPPA to the host data, which here is assumed to have acceleration in
% the first phase-encoding direction only (as is the current limitation of
% the MP2RAGE Siemens code).
%
% This version attempts to loop over 'Echoes' and 'Sets'. Note that if there
% are multiple 'Repetitions' then only the repetition specified by
% 'iRepToRecon' will be reconstructed - and it is assumed that there is
% separate ACS data available for each repetition.
%
%
% -------
% daniel.gallichan@epfl.ch, June 2015
%
% 11/3/16 - danielg - created '_RAMonly' version which should handle more
%                     arbitrary datasets, as when using matfiles it is
%                     difficult to handle arrays which might have a different
%                     number of dimensions. Also removed the coil
%                     compression as it doesn't seem to be compatible with
%                     motion-correction. As it's RAM only the mOutGRAPPA
%                     structure is also only constructed at the end, as then 
%                     it becomes compatible with parfor...
   

if nargin < 4
    nSliceNeighbours = 2; % no. of neighbouring 'virtual' slices to include in ACS lines to try to improve conditioning and regularize in the readout direction
    % quick tests suggest that 2 is a good number to use... :)
end
if nargin < 2
    iRepToRecon = 1;
end

if nargin < 3
    iEcoToRecon = 1;
end

%%%%%%%%%%%%%%%%%%%%
%%% Will only currently work with 1D GRAPPA...!
%%%%%%%%%%%%%%%%%%%%
%%


figIndex = 999; % use this figure number for live 'show' of how things are going...

nread = twix_obj.image.dataSize(1);
npe1 = twix_obj.image.dataSize(3);
nsli = twix_obj.image.dataSize(5);
nc = twix_obj.image.dataSize(2);

nSet = twix_obj.image.dataSize(10);
nEco = twix_obj.image.NEco;

nRep = twix_obj.image.NRep;

sliceIndices = twix_obj.hdr.Config.chronSliceIndices(1:nSli);

[~, sliceIndicesSorted] = sort(sliceIndices, 'ascend');

                        
if iRepToRecon > nRep
    disp(['Error, tried to recon repetition ' num2str(iRepToRecon) ', but only found ' num2str(nRep) ' repetitions in data']);
    return
end

if iEcoToRecon > nEco
    disp(['Error, tried to recon repetition ' num2str(iEcoToRecon) ', but only found ' num2str(nEco) ' echoes in data']);
    return
end

nACSmeas1 = twix_obj.refscan.dataSize(3);

nACSmax = 40;
nACS1 = min(nACSmeas1,nACSmax);
nACS2 = 1;

centerLin1 = twix_obj.image.centerLin(1);

tic
disp('..............')
disp('... Loading the raw data file into RAM')
allData = twix_obj.image(:,:,:,:,:,:,:,iEcoToRecon,iRepToRecon,:);
allData = permute(allData, [5, 2, 3, 1, 4]);
allData = allData(sliceIndicesSorted,:,:,:,:,:);
disp('Done')
disp('..............')


%% Now do GRAPPA on each of those slices


ACSdata = twix_obj.refscan(:,:,:,:,:,:,:,iEcoToRecon,iRepToRecon,:);
ACSdata = permute(ACSdata, [5, 2, 3, 1, 4]);
ACSdata = ACSdata(sliceIndicesSorted,:,:,:,:,:);
%%% Define GRAPPA kernel
gx = 4; gy = 3; % size of GRAPPA kernel
Ax = twix_obj.hdr.MeasYaps.sPat.lAccelFactPE;
Ay = 1;
grapKernel = zeros(gx + (gx-1)*(Ax-1), gy + (gy-1) * (Ay-1));

% source points are marked with 1
grapKernel(1:Ax:end, 1:Ay:end) = 1;

% target points are marked with 0.5
if gx == 1,    startx = 2; else  startx = 2 + (floor(gx/2) - 1) * Ax; end
if gy == 1,    starty = 2; else  starty = 2 + (floor(gy/2) - 1) * Ay; end

grapKernel(startx:startx+Ax-2, starty:starty+Ay-2) = 0.5;
grapKernel(startx:startx+Ax-1, starty:starty+Ay-2) = 0.5;
grapKernel(startx:startx+Ax-2, starty:starty+Ay-1) = 0.5;
%%%


disp('...............')
disp('Declaring temporary variables to store GRAPPA recon')
tic
grappaRecon = complex(zeros(nsli,nc,npe1,nread,nSet,'single'));
% reconSoS = zeros(nread,npe1,npe2,nSet,'single');

timingReport.declareVariables = toc;
disp('Done')
disp('...............')

    
%%
tic
disp('..............')
disp('...Now doing GRAPPA on each ''slice'' in the readout direction...')

parfor iS = 1:nsli % was nread
        
        for iSet = 1:nSet

            iPE2 = [1:nACS2]-floor(nACS2/2)+centerPar1;
            if any(iPE2<1 | iPE2>nACSmeas2) % seems this can happen for certain orientations / choices of parameters...
                iPE2 = 1:nACS2; % then hope that this works in this case...
            end
                        
            thisdata = squeeze(allData(iS,:,:,:,1,1,1,1,1,iSet));   % already chosen iEco        
            thisdata = permute(thisdata,[2 3 1]);
            
            startX = find(thisdata(:,1,1),1,'first');            
           
            if nSliceNeighbours > 0
                
                src = []; targ = [];
                
                for this_iS = iS-nSliceNeighbours:iS+nSliceNeighbours
                    if this_iS > 1 && this_iS <= nread
                        thisACSdataFull = squeeze(ACSdata(this_iS,:,:,:,1,1,1,1,1,iSet)); % iEcoToRecon already set
                        thisACSdataFull = permute(thisACSdataFull,[2 3 1]);
                        thisACSdata = thisACSdataFull([1:nACS1]-floor(nACS1/2)+floor(nACSmeas1/2),iPE2,:);
                        [~, this_src, this_targ] = GrappaCalib3D_arb(thisACSdata,grapKernel,0,0);
                        src = [src; this_src];
                        targ = [targ; this_targ];
                    end
                end
                
                % and store the current slice for reinsertion
                thisACSdataFull = squeeze(ACSdata(iS,:,:,:,1,1,1,1,1,iSet)); % iEcoToRecon already set
                thisACSdataFull = permute(thisACSdataFull,[2 3 1]);
                
                thisGrapW = pinv(src)*targ;
                          
            else
            
                % no slice neighbors
                thisACSdataFull = squeeze(ACSdata(iS,:,:,:,1,1,1,1,1,iSet)); % iEcoToRecon already set
                thisACSdataFull = permute(thisACSdataFull,[2 3 1]); % this permute is necessary for my GRAPPA code, but may not be the fastest approach...
                thisACSdata = thisACSdataFull([1:nACS1]-floor(nACS1/2)+floor(nACSmeas1/2),iPE2,:);

                thisGrapW = GrappaCalib3D_arb(thisACSdata,grapKernel,0,1);
        
            end
            
            thisGrapRecon = GrappaReco3D_arb(thisdata,grapKernel,thisGrapW,[Ax 1],[startX 1]);
            
            
            %
            % and put the ACS lines back in
            % WARNING - I don't know if the following line will be correct for all
            % possible image resolution/matrix sizes, so possibly safer to leave this to the
            % separate image recon part...
            thisGrapRecon([1:nACSmeas1]-floor(nACSmeas1/2)+centerLin1-1,:,:) = thisACSdataFull;
            
                        
            grappaRecon(iS,:,:,:,iSet) = reshape(permute(thisGrapRecon,[3 1 2]),[1 nc npe1 nread]);

            
%             if ~isempty(figIndex)
%                 fig(figIndex)
%                 clf
%                 imab(thisImage)
%                 title(['GRAPPA recon, slice ' num2str(iS) ' out of ' num2str(nread) ', echo index: ' num2str(iEco) ', set index: ' num2str(iSet)])
%                 drawnow
%             end
            
            
        end % iSet loop
        
    %end %iEco loop
    
    if mod(iS,10)==0
        fprintf('.');
    end
    
end % iS for 'slices' in readout direction

fprintf('\n');
disp('Done')
disp('..............')
timingReport.GRAPPArecon = toc;

mOutGRAPPA.reconSoS = reconSoS;
mOutGRAPPA.timingReport = timingReport;
if ~isempty(combinePars)
    mOutGRAPPA.dataCombined = dataCombined;
end


