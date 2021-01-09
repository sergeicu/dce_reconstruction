% ## Reshape the given 1D vector into a complex 2D image. The real part
% ## is taken from the first half of the vector and the imaginary part from
% ## the second half.

function [img]=vec_to_img(vec,br)
    img=reshape(vec(1:end/2),[br, br]) + reshape(1i*vec(end/2+1:end),[br, br]);
end
