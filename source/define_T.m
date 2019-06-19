function [Ei,Et,t,gamma,thresh,Nx,Ny,NLx,NLy,OffsetX,OffsetY,steps] = define_T()
% [T,Et,circ_in,circ_outside,domain] = define_T()
% Function to define T = IFFT[Et*t]
%
% input = none
%
% Outputs:
% Ei, incidente complex input field
% Et, target field
% t, target domain where Et is active
% gamma, exterior set of t
% thresh, value to stop the phase retrieval process
% Nx,Ny, dimension of the calculated mask
% NLx,NLy, dimension of the target field
% OffsetX,OffsetY, displacement of the origin for order zero problem
% steps, maximum number of iteration
%
% Date: 05/30/2019
% Authors: A. Federico - M. Yommi
%

%% Input conditions
    % Maximum mplitude of the Incident field
        Amp_Ei = 24;  % alpha value relationship

    % Size of the target field Er
        NLy = 128;
        NLx = NLy;

    % Size of the phase-only mask tau  
        Nx = 270;
        Ny = 320;

    % Target field
        OffsetX = 0;  % 72; % zero order compensation
        OffsetY = 0;  % 97;
        
    % Incident random field generation (speckle field)
        rng(33); % set seed
        U = speckle_gen(135,Ny,Nx);
        as = max(abs(U(:)));
        U = U/as;
        Ei = Amp_Ei.*U;    
        
    % Domains and field target definitions
        uno = ones(Ny,Nx); 
        t=zeros(Ny,Nx);
        t(Ny/2-NLy/2+OffsetY:Ny/2+NLy/2-1+OffsetY,Nx/2-NLx/2+OffsetX:Nx/2+NLx/2-1+OffsetX) = 1;
        gamma = uno-t;

        A_t = zeros(Ny,Nx);
        L = imread('../img/baboon.tif','tif');
        L = imresize(L,[NLy NLx]); % [numrows numcols]
        baboon = mat2gray(double(L));

        Phi_t = zeros(Ny,Nx);
        L = imread('../img/lenna.tif','tif');
        L = imresize(L,[NLy NLx]);  % [numrows numcols])
        Lenna = mat2gray(double(L))*2*pi-pi; 

    % Target field
        YY = baboon.*exp(1i*Lenna);

        A_t(Ny/2-NLy/2+OffsetY:Ny/2+NLy/2-1+OffsetY,Nx/2-NLx/2+OffsetX:Nx/2+NLx/2-1+OffsetX) = abs(YY);
        Phi_t(Ny/2-NLy/2+OffsetY:Ny/2+NLy/2-1+OffsetY,Nx/2-NLx/2+OffsetX:Nx/2+NLx/2-1+OffsetX) = angle(YY);
    
        Et = A_t.*exp(1i*Phi_t); % Target beam embedded in Ny by Nx pixels
        
    % Threshold settings
        thresh = 2e-6; % phase retrieval comparison
        
    % Maximum iteration number
        steps=1000;
        
return
