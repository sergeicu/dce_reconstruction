%##########################################################################
%
%  Matlab source code for MRI reconstruction method described in:
%  
%  Block KT, Uecker M, Frahm J. 
%  Undersampled radial MRI with multiple coils. Iterative image 
%  reconstruction using a total variation constraint. 
%  Magn Reson Med. 2007 Jun;57(6):1086-98
% 
%  The algorithm in this file corresponds to Figures 6-8.
%
%  The source code uses the following external packages:
%    - NUFFT toolkit by Jeffrey Fessler 
%      (http://www.eecs.umich.edu/~fessler/)
%    - NUFFT operator by Miki Lustig 
%      (http://www.eecs.berkeley.edu/~mlustig/Software.html)
%    - Siemens TWIX file reader by Philipp Ehses
%    - Poblano Toolbox by Sandia National Laboratories 
%      (http://software.sandia.gov/trac/poblano)
%    - MRI Phantom by Ronald Ouwekerk
%      (http://www.mathworks.com/matlabcentral/fileexchange/1759-mriphantom)
%
%  If you are using this code, please cite the publication listed above.
%
%  Version 25.03.15.
%
%##########################################################################

% ## Clear previous data and variables
% clc; 
% close all;
% clear all;
function [k, coilprofile, dcf,gridding]= calc_iterativeradial_multicoil(datax,chunks);
    % Set to 0 for phantom data and 1 for brain data
    dataBrain=1;
    display('vesrion')
    % datax=datax;
    % ## Set reconstruction parameters

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

    %param.filename='data/phantom.dat';

    if (dataBrain==1)
    %   param.filename='data/brain.dat';
        param.lambdaImage=0.000005; % for brain
        param.lambdaCoils=15;
        param.stopTolImage=1e-8;
        gradCorr=0.365;
    end


    % ## Include packages in subfolders
    % addpath('../nufft');                % for the NUFFT
    % addpath('../utilities');            % for the NUFFT
    % addpath('../graph');                % for the NUFFT
    % addpath('../imagescn_R2008a/');     % for plotting images
    % addpath('../read_meas_VBVD/');      % for reading TWIX files
    % addpath('../poblano_toolbox_1.1');  % for numerical optimization
    % addpath('../mriphantom');           % for phantom simulation


    % ## Read the data from the Siemens TWIX file
    %data=mapVBVD(param.filename);
    %load sila;
    % data=double(datax);

    % For VD software, the TWIX file may contain adjustmentdata. We only want
    % to look at the image data for now.
    %data=data.image;

    % % Read the sequence params from the file
    % baseresolution=double(data.NCol);
    % spokes=double(data.NLin);
    % channels=double(data.NCha);
    % slices=1;
    [baseresolution channels spokes]=size(datax);
    %keyboard
    % baseresolution=448;%double(data.NCol);
    % spokes=930;%double(data.NLin);
    % channels=34;%double(data.NCha);
    slices=1;


    fprintf('\n');
    fprintf('Radial Iterative Reconstruction\n');
    fprintf('-------------------------------\n\n');
    fprintf('Spokes = %d\n',         spokes        );
    fprintf('Channels = %d\n',       channels      );
    fprintf('Baseresolution = %d\n', baseresolution);
    fprintf('\n');

    % Get the k-space data. Data comes in as [samples,channels,spokes,slice]
    %rawdata=data().image;
    rawdata=double(datax);
    clear datax;

    % ## Downsampling of readout data (reduces matrix sizes)
    if (param.readoutDownsampling==1)
        tempraw=fftshift(fft(ifftshift(rawdata(:,:,:,1),1),[],1),1);
        tempraw=tempraw(end/4+1:3*end/4,:,:,:);
        rawdata=ifftshift(ifft(fftshift(tempraw,1),[],1),1);
        baseresolution=baseresolution/2;
    end

    spokesOffset=0;

    % ## Overwrite sequence params for retrospective undersampling
    if (param.undersampleSpokesTo~=0)
        % Rejecting spokes at the beginning rather than at then end to reduce
        % steady-state effects that affect the initial spokes
        spokesOffset=spokes-param.undersampleSpokesTo;
        spokes=double(param.undersampleSpokesTo);
        fprintf('Undersampling spokes to = %d\n', spokes);
    end


    % ## Prepare the density compensation function (DCF)
    dcfRow=ones(baseresolution,1);
    for i=1:baseresolution
    dcfRow(i)=abs(baseresolution/2-(i-0.5));
    end
    dcfRow=pi/spokes*dcfRow;
    dcf=repmat(dcfRow,1,spokes);


    dcf1=repmat(dcfRow/chunks,1,spokes/chunks);

    clear dcfRow;

    spokeschunk1=spokes/chunks;
    % ## Calculate angles for the radial trajectory
    if (param.orderingGoldenAngle==1)
        % For the Golden-Angle scheme
        GA = 111.246117975/180*pi; 
        phi = [pi/2+spokesOffset*GA:GA:GA*spokeschunk1+spokesOffset*GA];
    else
        % For the linear ordering mode
        phi = zeros(1,spokeschunk1);
        for i=1:spokeschunk1
            phi(i)=pi/spokeschunk1*(i-1);
            % Flip every second spoke, as done in the sequence
            if (mod(i,2) == 0)
                phi(i)=phi(i)+pi;
            end
        end
    end
    phi = mod(phi,2*pi);
    phi1=[];
    for ii=1:chunks
        phi1=[phi1 phi];
    end
    phi=phi1;

    clear phi1;
    % ## Calculate the k-space trajectory points 
    if (param.readoutDownsampling==0)
        % Without readout downsampling
        rho = linspace(0,baseresolution-1,baseresolution)' - (baseresolution-1)/2;
    else
        % With readout downsampling
        % gradCorr contains a dataset dependend shift to correct for gradient
        % delays (would be better to estimate them from the data)
        rho = linspace(0,baseresolution-1,baseresolution)' - (baseresolution-gradCorr)/2;
    end

    % Norm to range -0.5 to 0.5 for the NUFFT
    rho = rho/baseresolution;

    % Generate vector with k-space coordinates (as complex values kx + i ky)
    k = double(rho*exp(-1j*phi));
    k1=double(rho*exp(-1j*phi(1:length(phi)/chunks)));

    % ## Simulate k-space data with numerical phantom (will overwrite
    % measurement data with numerical Shepp-Logan data)
    if (param.simulateData==1)
        fprintf('Simulating data...\n');
        phantomdata=mriphantom(imag(-1*k*baseresolution/(0.7*pi)),real(-1*k*baseresolution/(0.7*pi)));
        for ic=1:channels
            rawdata(:,ic,1:spokes,1)=phantomdata;
        end
        fprintf('Done.\n');
    end


    % ## Prepare the NUFFT
    fftSize=[baseresolution,baseresolution];
    FT = NUFFT(k1, 1, 1, [0,0], fftSize, 2);


    % ## Perform gridding reconstruction for comparison
    fprintf('Calculating gridding reconstruction...');

    % Create image to sum the gridded channels (for sum-of-squares combination)
    gridding=zeros(baseresolution,baseresolution);
    % display('')   
    % Loop over channels and calculate image for each coil
    for channel=1:channels
        % Fetch rawdata for slice and channel and multiply with density
        % compensation function
        workRawData=dcf1.*double(squeeze(rawdata(:,channel,1+spokesOffset:spokes/chunks+spokesOffset,1)));

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
    param.rawdata=rawdata(:,:,1:spokes/chunks);
    param.br=baseresolution;
    param.spokes=spokes/chunks;
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
    options.Display='off';

    % keyboard
    % Loop over the channels to create "smooth" image for each coil
    for ic=1:channels
        % Read k-space data for the channel 
        param.y=double(squeeze(rawdata(:,ic,1+spokesOffset:spokes/chunks+spokesOffset,1)));

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


    % % ## Step 2: Calculation of final joint-coil image
    % 
    % % Initialize the optimizer
    % param.coilprofile=coilprofile;
    % param.y=double(squeeze(rawdata(:,:,1+spokesOffset:spokes+spokesOffset,1)));
    % x0=zeros(param.nvar,1);
    % iterationCounter=0;
    % options.StopTol=param.stopTolImage;
    % 
    % % First run the optimizer for 10 iterations without penalty terms
    % param.enablePenalty=0;
    % options.MaxIters=10;
    % out = lbfgs(@(x) costfunction_image(x,param), x0, options);
    % 
    % % Now enable the penalty terms and run the optimizer for the remaining
    % % iterations.
    % param.enablePenalty=1;
    % options.MaxIters=param.iterationsImage-10;
    % out = lbfgs(@(x) costfunction_image(x,param), out.X, options);
    % 
    % % Reshape the result vector into 2D image format
    % iterative=vec_to_img(out.X,param.br);


    % % ## Now show the results
    % fprintf('Finished.');
    % close;
    % 
    % % Get magnitude and crop center part of images. Show the reconstructed 
    % % images with 4x bilinear interpolation (looks better)
    % figure('Name', 'Coil Profiles'),     imagescn(imresize(fliplr(abs(coilprofile)),[size(coilprofile,1)*4 size(coilprofile,2)*4],'bilinear'),[],[],[],3);
    % figure('Name', 'Gridding Solution'), imagescn(imresize(fliplr(abs(gridding)),   [size(gridding,1)*4 size(gridding,2)*4],      'bilinear'),[],[],[],3);
    % figure('Name', 'Iterative Solution'),imagescn(imresize(fliplr(abs(iterative)),  [size(iterative,1)*4 size(iterative,2)*4],    'bilinear'),[],[],[],3);
end