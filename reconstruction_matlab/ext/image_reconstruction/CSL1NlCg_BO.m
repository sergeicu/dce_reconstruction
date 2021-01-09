function x = CSL1NlCg_BO(x0,param)
% 
% res = CSL1NlCg(param)
%
% Compressed sensing reconstruction of undersampled k-space MRI data
%
% L1-norm minimization using non linear conjugate gradient iterations
% 
% Given the acquisition model y = E*x, and the sparsifying transform W, 
% the pogram finds the x that minimizes the following objective function:
%
% f(x) = ||E*x - y||^2 + lambda * ||W*x||_1 
%
% Based on the paper: Sparse MRI: The application of compressed sensing for rapid MR imaging. 
% Lustig M, Donoho D, Pauly JM. Magn Reson Med. 2007 Dec;58(6):1182-95.
%
% Ricardo Otazo, NYU 2008
%

    fprintf('\n Non-linear conjugate gradient algorithm')
    fprintf('\n ---------------------------------------------\n')

    % starting point
    x=x0;

    % line search parameters
    maxlsiter = 150 ;
    gradToll = 1e-3 ;
    param.l1Smooth = 1e-15;	
    alpha = 0.01;  
    beta = 0.5;
    t0 = 1 ; 
    k = 0;

    normalization = 1;
       
    GP_params = struct();

    GP_params.sigma_f = 1;
    GP_params.sigma_n = 1e-6;
    GP_params.lb = 0;
    GP_params.ub = 1;
    GP_params.length = GP_params.ub/4;
    

    % compute g0  = grad(f(x))
    g0 = grad(x,param);
    dx = -g0;
    
    normY = norm(param.y(:))^2;

    % iterations
    while(1)

        % backtracking line-search
        fcn = @(t) -objective(x,dx,t,param);
        if k == 0
            f1 = fcn(0);
        end
        t0 = [0, GP_params.ub];
        f0 = [f1, fcn(t0(end))];
        [t, f1] = bayesOpt_GaussianProcess(t0,f0,fcn,GP_params,'');

        % control the number of line searches by adapting the initial step search
        if abs( t - GP_params.ub) < 1e-6
            GP_params.ub = GP_params.ub / beta; 
            GP_params.length = GP_params.ub/4;
        else
            GP_params.ub = GP_params.ub * beta;
            GP_params.length = GP_params.ub/4;
        end

        % update x
        x = (x + t*dx);

        %conjugate gradient calculation
        g1 = grad(x,param);
        bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
        g0 = g1;
        dx =  - g1 + bk* dx;
        k = k + 1;
        
        % print some numbers	
        if param.display,
            fprintf(' ite = %d, cost = %0.6f, grad = %0.6f \n',k, -f1/normY, norm(dx(:))/normY );
        end

        % stopping criteria (to be improved)
BenvinGuda        if (k > param.nite) || (norm(dx(:)) < gradToll), break;end
        
        normalization = abs(f1);
%         close all;
    end
    return;
end
function res = objective(x,dx,t,param) %**********************************

    % L2-norm part
    w=param.E*(x+t*dx)-param.y;
    L2Obj=w(:)'*w(:);

    % L1-norm part
    if param.lambda
       w = param.W*(x+t*dx); 
       L1Obj = sum((conj(w(:)).*w(:)+param.l1Smooth).^(1/2));
    else
        L1Obj=0;
    end

    % objective function
    res = L2Obj+param.lambda*L1Obj;
end
function g = grad(x,param)%***********************************************

    % L2-norm part
    L2Grad = 2.*(param.E'*(param.E*x-param.y));

    % L1-norm part
    if param.lambda
        w = param.W*x;
        L1Grad = param.W'*(w.*(w.*conj(w)+param.l1Smooth).^(-0.5));
    else
        L1Grad=0;
    end


    % composite gradient
    g=L2Grad+param.lambda*L1Grad;

end

