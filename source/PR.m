function [tau,rms_angle,rms_fourier] = PR(Ei,Et,thresh,gamma,t,NLx,NLy,steps)
% 
% Find tau for input arguments: Ei,Et,thresh,gamma,t,NLx,NLy
%
% input = Ei,Et,thresh,gamma,t,NLx,NLy
%
% outputs = tau (phase mask), RMS values: rms_angle (phase),rms_fourier(intensity)
% 
%
% Date: 10/30/2019
% Authors: A. Federico - M. Yommi
%
    T = iFFT2(Et);
    rng(222); phi = rand(size(T))*2*pi-pi; 
    rms_angle(1) = 1; 
    k = 0; 
    Ai = abs(Ei);
    
    h = waitbar(0,'Please wait...');
    while rms_angle > thresh
        k=k+1;  % number of iterations
        if k == steps 
            break
        end

        waitbar(k / steps);
                        
        c = T.*exp(-1i*phi)+conj(T).*exp(1i*phi);
        B0 = (-c+sqrt(c.^2-4*(abs(T).^2-Ai.^2)))/2;  % positive solution

        B = iFFT2(FFT2(T+B0.*exp(1i*phi)).*gamma);
        phi = angle(B);
        tau = (T+B0.*exp(1i*phi))./Ei;

        FTt = FFT2(tau.*Ei);
        S = FTt.*t;
        rms_fourier(k) = sum((abs(Et(:))-abs(S(:))).^2)/NLy/NLx;         
        rms_angle(k) = sum((angle(S(:))-angle(Et(:))).^2)/NLy/NLx;
    end
    
    close(h)
    
    save('../env/B.mat', 'B');
    save('../env/tau.mat', 'tau');
end