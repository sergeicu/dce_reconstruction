function  res = TV_multiscale(x_prev)

    res.adjoint = 0;
    res.x_prev = x_prev;
    res.T = size(x_prev,3);
    res = class(res,'TV_multiscale');

end