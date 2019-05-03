% Homework 4
% Delay-and-sum Beamformer, using weights

    close all
    clear
    clc

    addpath('..\tools')
    
%% Parameters:
    theta   = 20;                   % direction of arrival (degrees)
    doa     = theta * pi/180;       % vector of directions of arrival (radians)
    winSel  = 9;
    
    f      = 1.0e3;                 % carrier frequency (narrow-band signal)

    c       = soundSpeed();         % speed of sound (m/s, t = 20°C by default)
    lambda  = c./f;                 % wavelength (meters)

    M       = 16;                   % number of sensors
    D       = lambda/2;             % separation of consecutive sensors (meters)
                                    % for a standard ULA, the inter-element spacing is half a wavelength
                                
    d       = D/lambda;             % this ratio should be less than 0.5 (grating lobes may appear otherwise)
    
%% Uniform linear array(ULA):    
    switch winSel
        case 1
            fhandle = @rectwin;
            disp('Rectangular window')
        case 2
            fhandle = @bartlett;
            disp('Triangular window')
        case 3
            fhandle = @hann;
            disp('Hann window')
        case 4
            fhandle = @hamming;
            disp('Hamming window')
        case 5
            fhandle = @blackman;
            disp('Blackman window')
        case 6
            fhandle = @blackmanharris;
            disp('Blackman-Harris window')
        case 7
            fhandle = @chebwin;
            disp('Chebyshev window')
        case 8
            fhandle = @gausswin;
            disp('Gaussian window')
        case 9
            fhandle = @taylorwin;
            disp('Taylor window')
    end
            
    
    w = window(fhandle, M);
    
%% Beam pattern:
    N_theta     = 1000;
    doa_theta   = linspace(-pi/2, pi/2, N_theta);       % general DOA vector 
    G_theta     = zeros(N_theta, 1);
    m           = (0:M-1)';                             % indices of sensors
    for th = 1:N_theta
        aULA_theta = exp(2j*pi*d*sin(doa_theta(th) - doa)*m); % general ULA vector
        G_theta(th) = w'*aULA_theta;
    end    
    
%% Display figures:
    
    figure('units','normalized','outerposition',[0 0 1 1])
    subplot 121
    plot(doa_theta*180/pi, 10*log10(abs(G_theta)/max(abs(G_theta))), 'LineWidth', 2)
    axis([-90 90 -60 0])
    xlabel 'DOA (degrees)'
    title 'Beam pattern'

    subplot 122
    polarplot(doa_theta, abs(G_theta)/max(abs(G_theta)), 'LineWidth', 2)
    ax = gca;
    ax.ThetaLim = [-90 90];