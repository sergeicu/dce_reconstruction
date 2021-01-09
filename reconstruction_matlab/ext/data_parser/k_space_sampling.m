%%
%
%       This function returns the position of the K-space samples from given protocols.
%
%
%       Available protocols:
%
%           - Stack of stars: 'stack_stars' - 
%                   - 
%
%
%       INPUT: 
%           - samplingType - str - type of sampling used for K-space.
%           - N     - int -  number of samples per spoke.
%           - NCha  - int -  number of chanels.
%           - NSpk  - int -  number of spokes.
%           - NChunk -int -  number of blocks with which the data was acquired.
%
%
%
%       OUTPUT:
%           - k_samples - <baseresolution,1>double -    k space coordinates ascomplex values (kx + i ky)
%           - dcf       - <baseresolution,Nspk,NSli>double - density compensation function
%
%   
%       AUHTOR:
%           Jaume Coll-Font <jcollfont@gmail.com>  (wrapping of Sila Kurugol's code)
%       
%
function [k_samples,dcf] = k_space_sampling( samplingType, N, NCha, NSpk, NSli, NChunk )

    if strcmp(samplingType , 'stack_stars' )
        [k_samples,dcf]= calc_k_dcf_stack_of_stars( N, NCha, NSpk, NChunk);
        dcf = repmat(dcf,[1,1,NSli]);
        k_samples = repmat(k_samples,[1,1,NSli]);
    end
end




function [k,dcf]= calc_k_dcf_stack_of_stars( baseresolution, channels, spokes, chunks)
    %for all chunks
    spokeschunk1=spokes/chunks;
    slices=1;
    param.orderingGoldenAngle=1;
    
    fprintf('\n');
    fprintf('Radial Iterative Reconstruction\n');
    fprintf('-------------------------------\n\n');
    fprintf('Spokes = %d\n',         spokes        );
    fprintf('Channels = %d\n',       channels      );
    fprintf('Baseresolution = %d\n', baseresolution);
    fprintf('\n');
    
    
    
    spokesOffset=0;
    gradCorr=0.4;
    
    % ## Prepare the density compensation function (DCF)
    dcfRow=ones(baseresolution,1);
    for i=1:baseresolution
       dcfRow(i)=abs(baseresolution/2-(i-0.5));
    end
    dcfRow=pi/spokes*dcfRow;
    dcf=repmat(dcfRow,1,spokes);
    clear dcfRow;
    
    
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
    
    rho = linspace(0,baseresolution-1,baseresolution)' - (baseresolution-gradCorr)/2;
    
    
    % Norm to range -0.5 to 0.5 for the NUFFT
    rho = rho/baseresolution;
    
    % Generate vector with k-space coordinates (as complex values kx + i ky)
    k = double(rho*exp(-1j*phi));
    
end