function coilprofile= calc_iterativeradial_multicoil_onechunk(datax)
% Note: The penalty weights lambdaCoils / lambdaImage and the stop criterion
% for the iterations stopTolCoils / stopTolImage need to be adjusted
% depending on the value range of the data, and on the undersampling level
param.lambdaCoils=10;
param.stopTolCoils=1e-7;
param.iterationsCoils=20;% was 60

param.lambdaImage=0.00002; % for phantom
param.iterationsImage=60;
param.stopTolImage=4e-8;

% Retrospective undersampling of golden-angle data (0=full data set)
param.undersampleSpokesTo=0;
% Ordering scheme (1=Golden Angle, 0=linear)
param.orderingGoldenAngle=1;
% Downsampling of readout oversampling to reduce matrix size (1=on, 0=off)
param.readoutDownsampling=0;
% Numerical simulation of k-space data (will overwrite scan data)
param.simulateData=0;
% Gradient-delay correction factor
gradCorr=0.4;


[baseresolution,channels,spokes]=size(datax);
%rawdata=data().image;


rawdata=double(datax);
clear datax;



% ## Prepare the density compensation function (DCF)
dcfRow=ones(baseresolution,1);
for i=1:baseresolution
   dcfRow(i)=abs(baseresolution/2-(i-0.5));
end
dcfRow=pi/spokes*dcfRow;
dcf1=repmat(dcfRow,1,spokes);
clear dcfRow;

spokesOffset=0;

if (param.orderingGoldenAngle==1)
    % For the Golden-Angle scheme
    GA = 111.246117975/180*pi; 
    phi = [pi/2+spokesOffset*GA:GA:GA*spokes+spokesOffset*GA];
else
    % For the linear ordering mode
    phi = zeros(1,spokes);
    for i=1:spokeschunk1
        phi(i)=pi/spokes*(i-1);
        % Flip every second spoke, as done in the sequence
        if (mod(i,2) == 0)
            phi(i)=phi(i)+pi;
        end
    end
end
phi = mod(phi,2*pi);

rho = linspace(0,baseresolution-1,baseresolution)' - (baseresolution-1)/2;


% Norm to range -0.5 to 0.5 for the NUFFT
rho = rho/baseresolution;

% Generate vector with k-space coordinates (as complex values kx + i ky)
k = double(rho*exp(-1j*phi));


% ## Prepare the NUFFT
fftSize=[baseresolution,baseresolution];
FT = NUFFT(k, 1, 1, [0,0], fftSize, 2);


% ## Perform gridding reconstruction for comparison
fprintf('Calculating gridding reconstruction...');

% Create image to sum the gridded channels (for sum-of-squares combination)
gridding=zeros(baseresolution,baseresolution);
% display('')   
% Loop over channels and calculate image for each coil
for channel=1:channels
    % Fetch rawdata for slice and channel and multiply with density
    % compensation function
    workRawData=dcf1.*double(squeeze(rawdata(:,channel,1+spokesOffset:spokes+spokesOffset,1)));

    % Run the NUFFT
    workImage=FT'*(workRawData);

    % Add squared channel image to buffer
    gridding = gridding + abs(workImage.*workImage);
end

% Calculate the root (as final part of the sum-of-squares calculation)
gridding = sqrt(gridding);
fprintf('done.\n');
% 
% 
% ## Iterative reconstruction

% Pass information to the optimizer / cost function
param.nvar=2*baseresolution^2;
param.FT=FT;
param.dcf=dcf1;
param.rawdata=rawdata;
param.br=baseresolution;
param.spokes=spokes;
param.channels=channels;
clearvars FT dcf1;
% Global counter for displaying the number of evaluations
global iterationCounter
iterationCounter=0;

% ## Step 1: Estimation of coil profiles

% Initialize array for storing the coil profiles
coilprofile=zeros(baseresolution,baseresolution,channels);

% Create figure for displaying the iterations
% figure('Name', 'Iterations')
options.Display='final';
% Loop over the channels to create "smooth" image for each coil
for ic=1:channels
    % Read k-space data for the channel 
    param.y=double(squeeze(rawdata(:,ic,1+spokesOffset:spokes+spokesOffset,1)));

    % Initialize optimizer with empty image
    x0=zeros(param.nvar,1);
    iterationCounter=0;

%     param.rawdata=rawdata;
    % Set stop criterion depending on the value range of the raw data
    options.StopTol=param.stopTolCoils;

    % Run the optimizer for 10 iteration without penalty terms
    param.enablePenalty=0;
    options.MaxIters=10;
    out = lbfgs(@(x) costfunction_coils(x,param), x0, options);
    
    % Now enable the penalty terms and run the optimizer for the remaining
    % iterations
    param.enablePenalty=1;
    options.MaxIters=param.iterationsCoils-10;
    out = lbfgs(@(x) costfunction_coils(x,param), out.X, options);

%       keyboard;
    % Reshape the result vector into a 2D image and store it
    coilprofile(:,:,ic)=vec_to_img(out.X,param.br);
end
clear out;

clear param;

% Sum-of-squares calculation of the coil profiles
ssqcoil=zeros(baseresolution,baseresolution);
for ic=1:channels
    ssqcoil=ssqcoil+abs(coilprofile(:,:,ic).*coilprofile(:,:,ic));
end

ssqcoil=sqrt(ssqcoil);
for ic=1:channels
    coilprofile(:,:,ic)=coilprofile(:,:,ic)./ssqcoil;
end   
