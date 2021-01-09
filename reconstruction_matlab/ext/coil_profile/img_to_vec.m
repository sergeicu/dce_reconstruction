% ## Reshape the given complex 2D image into a 1D vector, with the real part
% ## being in the first half of the vector and the imaginary part in the 
% ## second half.

function [vec]=img_to_vec(img)
    vec=cat(1,reshape(real(img),[],1),reshape(imag(img),[],1));
end
