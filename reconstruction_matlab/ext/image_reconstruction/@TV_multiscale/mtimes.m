function res = mtimes(a,b)

    T = size(a.x_prev,3);
    ix = 1:T;%reshape(repmat(1:T,[2,1]),[1,2*T]);
    if a.adjoint
        res= b - a.x_prev(:,:,ix);
    else
        res = a.x_prev(:,:,ix)  - b;
        
    end
    
end