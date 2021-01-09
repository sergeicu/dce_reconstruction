% ## Helper function for calculating the L2 norm of the finite differences
% ## and its gradient for a given 2D complex image, which will smooth the 
% ## image. The penalty is calculated independeltly for the real and 
% ## imaginary part, using first-order finite differences.

function [f,g]=L2D(img)
    % Assuming that the image is squared
    br=size(img,1);

    % Initialize the output variables (f is double, and g complex matrix of
    % size br x br
    f=0;
    g=complex(zeros(br),0);
    
    % Loop over image and simultanously calculate the penalty and the
    % gradient. The edge pixels are excluded here for simplicity (a better
    % solution would be to "connect" the left-rigt and top-bottom edges)
    for ix=2:br-1
        for iy=2:br-1
            f=f+(real(img(ix,iy))-real(img(ix-1,iy)))^2;
            f=f+(real(img(ix,iy))-real(img(ix,iy-1)))^2;

            f=f+(imag(img(ix,iy))-imag(img(ix-1,iy)))^2;
            f=f+(imag(img(ix,iy))-imag(img(ix,iy-1)))^2;

            gradval=0;
            gradval=gradval+2*( real(img(ix,iy))  -real(img(ix-1,iy)) );
            gradval=gradval-2*( real(img(ix+1,iy))-real(img(ix,iy)) );

            gradval=gradval+2*( real(img(ix,iy))  -real(img(ix,iy-1)) );
            gradval=gradval-2*( real(img(ix,iy+1))-real(img(ix,iy)) );

            igradval=0;
            igradval=igradval+2*( imag(img(ix,iy))  -imag(img(ix-1,iy)) );
            igradval=igradval-2*( imag(img(ix+1,iy))-imag(img(ix,iy)) );

            igradval=igradval+2*( imag(img(ix,iy))  -imag(img(ix,iy-1)) );
            igradval=igradval-2*( imag(img(ix,iy+1))-imag(img(ix,iy)) );

            g(ix,iy)=gradval+ 1i*igradval;
        end
    end
    
end

