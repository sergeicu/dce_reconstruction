function res = mtimes(a,b)
    

    [Nx,Ny,T] = size(b);
    %% ASSUMING 2D images over time and fixed image sizes... everything is very well done!
    
    % reshape stuff
%     res = zeros(Nx*Ny,T);
%     
%     maskMask = zeros(Nx,Ny);
%     maskMask(Nx/4 + [1:Nx/2], Nx/4 + [1:Ny/2]) = 1;
%     maskMaskIX = find(maskMask(:));
    
    tempB = reshape(b,[Nx*Ny,T]);
    
    % compute differences
    res = tempB - a.LTIData;
    res = reshape(res,[Nx,Ny,T]);

end
