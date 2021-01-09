% ## Cost function for estimating the coil profiles.
%
% Input:
%   x     = current image estimate (in vector form)
%   param = reconstruction parameters
%
% Output:
%   f     = Cost function value
%   g     = Gradient of cost function

function [f,g]=costfunction_coils(x,param)

    global iterationCounter

    % Reshape the given estimate vector into 2D image form
    est_x=vec_to_img(x,param.br);
        
    % Show the current image estiamte
    %imagescn(real(est_x),[],[],[],3);
    iterationCounter=iterationCounter+1;
    %set(gcf, 'name', sprintf('Step 1 -- Coil Estimation Iterations -- (Evalutions: %d)',iterationCounter));
   % drawnow;  

    % ## Calculate the data-fidelity term (cost-function value and gradient)
    est_y=param.FT*est_x;
    res=est_y-param.y;  
    f=0.5*norm(res(:))^2;

    % Use preconditioning to accelerate convergence (needed for most
    % optimizers)
    res=param.dcf.*res;
    work=param.FT'*res;
    g=img_to_vec(work);

    % ## Calculate the penalty term
    if (param.enablePenalty==1)
        % Use helper function to calculate value and gradient of penalty
        [penaltyVal,penaltyGrad]=L2D(est_x);

        % Add to cost function
        %fprintf('Data fidelity=%f  Penalty=%f\n',f,param.lambdaCoils*penaltyVal);
        f=f + param.lambdaCoils*penaltyVal;
        g=g + param.lambdaCoils*img_to_vec(penaltyGrad);  
    end    
end

