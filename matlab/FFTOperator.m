function [A,At] = FFTOperator(mask)
% This function generates 3-D undersampling fourier operator
% DongWang
    [m,n,d] = size(mask);
    mask = fftshift(fftshift(mask,1),2);
%     mask = fftshift(mask);
    non_zero = find(mask~=0);
    A = @(x)FFTForward(x,non_zero);
    At = @(x)FFTBackward(x,non_zero);

    function A = FFTForward(img,non_zero)
        A = zeros(m,n,d);
        p = 1/sqrt(m*n)*fft2(img);
        A(non_zero) = p(non_zero);
%         A = A(:);
    end

    function At = FFTBackward(img,non_zero)
        p = zeros(m,n,d);
        p(non_zero) = img(non_zero);
        At = sqrt(m*n)*ifft2(p);
%         At = At(:);
    end

end
