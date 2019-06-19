function E = FFT2(E0)
% E = iFFT2(E0)
% Calculate the Fast Fourier Transform balanced
%
% input = E0
%
% output = E
%
% Date: 05/30/2019
% Authors: A. Federico - M. Yommi
%   

    [Ny,Nx] = size(E0);
    E = fftshift(fft2(ifftshift(E0)))/sqrt(Ny)/sqrt(Nx);
end    