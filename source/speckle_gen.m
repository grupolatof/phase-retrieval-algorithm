function U = speckle_gen(fsp,Ny,Nx)
%
% Function to generate an speckle field
%
% input = fsp, Ny, Nx. fsp is the radius of the frequency filter.
%
% output = random field U
%
% Date: 06/10/2019
% Authors: A. Federico - M. Yommi
%

    fi = rand(Ny,Nx)*2*pi-pi;
    U = exp(1i*fi);  
    IF = fftshift(fft2(U));
    [FU,FV] = freqspace(size(IF),'meshgrid');
    FU = FU*Nx/2;
    FV = FV*Ny/2;
    IF(sqrt(FU.^2+FV.^2)>fsp) = 0;  % lowpass 
    U = ifft2(ifftshift(IF));       % speckle field
end