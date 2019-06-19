function main ()
%% Initialization
    close all
           
    [Ei,Et,t,gamma,thresh,Nx,Ny,NLx,NLy,OffsetX,OffsetY,steps] = define_T();

%% Iterative algorithm for phase retrieval PR                   
    [tau,rms_angle,rms_fourier] = PR(Ei,Et,thresh,gamma,t,NLx,NLy,steps);
        
%% Figures of the incidente field
    figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
            s1 = subplot(1,2,1);
            imagesc(abs(Ei))
            axis square
            colorbar()
            format_subplot(s1)
            title('Amplitude of the incident light field Ei')
            s2 = subplot(1,2,2);
            imagesc(angle(Ei));
            title('Phase of the incident light field Ei')
            axis square
            colorbar()
            format_subplot(s2)
            
%% Results 8 bits discretization and comparison
% Discretization of the mask T
    alfa = 8; % Bits resolution of the computed mask
    fase_T_dis = double(uint8((2^alfa-1)/2/pi*angle(tau)+(2^alfa-1)/2))/(2^alfa-1);
    Er = FFT2(Ei.*exp(1i*(fase_T_dis*2*pi-pi)));

    figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
            s3 = subplot(1,2,1);
            imagesc(abs(Er),[0 1]);
            title('Amplitude of the recovered field Er')
            axis square
            colorbar()
            format_subplot(s3)
            s4 = subplot(1,2,2);
            imagesc(abs(Et));
            title('Amplitude of the target beam Et')
            axis square
            colorbar()
            format_subplot(s4)

    figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
            s5 = subplot(1,2,1);
            imagesc(angle(Er));
            title('Phase of the recovered field Er')
            axis square
            colorbar()
            format_subplot(s5)
            s6 = subplot(1,2,2);
            imagesc(angle(Et));
            title('Phase of the target beam Et')
            axis square
            colorbar()
            format_subplot(s6)
            
%% Figures of the determinated mask           
    figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
            s7 = subplot(1,2,1);
            imagesc(angle(tau));
            title('Phase of the phase mask Tau')
            axis square
            colorbar()
            format_subplot(s7)
            s8 = subplot(1,2,2);
            imagesc(abs(tau));
            title('Amplitude of the phase mask Tau')
            axis square
            colorbar()
            format_subplot(s8)
            
%% Convergence analysis
    figure('Color',[1 1 1],'units','normalized','outerposition',[0 0 1 1])
            s9 = subplot(1,2,1);
            plot(log10(rms_fourier))
            title('Evolution of the RMS values for the recovered amplitude')
            ylabel('rms_a')
            xlabel('Iteration number')
            format_subplot(s9)
            s10 = subplot(1,2,2);
            plot(log10(rms_angle))
            title('Evolution of the RMS values for the recovered phase')
            ylabel('rms_p')
            xlabel('Iteration number')
            format_subplot(s10)

%% Quality index determination
    FTt = FFT2(exp(1i*angle(tau)).*Ei); % focalization
                       
        alfa = 12; % CCD 12 bits discretization
        Intensity  = round(mat2gray(abs(Et(Ny/2-NLy/2+OffsetY:Ny/2+NLy/2-1+OffsetY,Nx/2-NLx/2+OffsetX:Nx/2+NLx/2-1+OffsetX)))*(2^alfa-1));
        Phase = round(mat2gray(angle(Et(Ny/2-NLy/2+OffsetY:Ny/2+NLy/2-1+OffsetY,Nx/2-NLx/2+OffsetX:Nx/2+NLx/2-1+OffsetX)))*(2^alfa-1));
        Intensity_recovered = round(mat2gray(abs(Er(Ny/2-NLy/2+OffsetY:Ny/2+NLy/2-1+OffsetY,Nx/2-NLx/2+OffsetX:Nx/2+NLx/2-1+OffsetX)))*(2^alfa-1));
        Phase_recovered = round(mat2gray(angle(Er(Ny/2-NLy/2+OffsetY:Ny/2+NLy/2-1+OffsetY,Nx/2-NLx/2+OffsetX:Nx/2+NLx/2-1+OffsetX)))*(2^alfa-1));

        fprintf('Discretization in the phase mask and CCD \n')
        
        [mssim,~] = ssim_index(Intensity , Intensity_recovered,[0.05 0.05],ones(5),2^alfa);
        fprintf('Quality_Intensity_for_phase_mask_8bits=%f8\n',mssim)
        [mssim,~] = ssim_index(Phase , Phase_recovered,[0.05 0.05],ones(5),2^alfa);
        fprintf('Quality_Phase_for_phase_mask_8bits=%f8\n',mssim)  
        
        Intensity_recovered = round(mat2gray(abs(FTt(Ny/2-NLy/2+OffsetY:Ny/2+NLy/2-1+OffsetY,Nx/2-NLx/2+OffsetX:Nx/2+NLx/2-1+OffsetX)))*(2^alfa-1));
        Phase_recovered = round(mat2gray(angle(FTt(Ny/2-NLy/2+OffsetY:Ny/2+NLy/2-1+OffsetY,Nx/2-NLx/2+OffsetX:Nx/2+NLx/2-1+OffsetX)))*(2^alfa-1));
        
        fprintf('\n')
        fprintf('Discretization in CCD\n')
        
        [mssim,~] = ssim_index(Intensity , Intensity_recovered,[0.05 0.05],ones(5),2^alfa-1);
        fprintf('Quality_Intensity=%f8\n',mssim)
        [mssim,~] = ssim_index(Phase , Phase_recovered,[0.05 0.05],ones(5),2^alfa-1);
        fprintf('Quality_Phase=%f8\n',mssim)
end
        
