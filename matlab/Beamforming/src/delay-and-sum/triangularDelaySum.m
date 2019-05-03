% Homework 5
% Triangular Delay-and-sum Beamformer varying its parameters (arbitrary triangle)

    close all
    clear
    clc

    addpath('..\tools')
    
%% Parameters:
    theta   = 20;                   % direction of arrival (degrees)
    doa     = theta * pi/180;       % vector of directions of arrival (radians)
    winSel  = 3;
    
    f      = 3.0e3;                 % carrier frequency (narrow-band signal)

    c       = soundSpeed();         % speed of sound (m/s, t = 20°C by default)
    lambda  = c/f;                  % wavelength (meters)

    M       = 3;                    % number of sensors (fixed to 3 for a triangle)
    D       = lambda/2;             % separation of consecutive sensors (meters)
                                    % for a standard ULA, the inter-element spacing is half a wavelength
                                
    d       = D/lambda;             % this ratio should be less than 0.5 (grating lobes may appear otherwise)
    
%% Reference system (assuming the reference sensor to be at: [x,y] = [0,0]):
    alpha = pi/3;           % angle where the third sensor is located (pi/3 for equilateral triangle)
    d21   = D;
    d31   = D;
    
    x1 = 0;                 y1 = 0;
    x2 = -d21;              y2 = 0;
    x3 = -d31*cos(alpha);   y3 = -d31*sin(alpha);    
    
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
    N_D         = 20;
    N_M         = 15;
    N_F         = 15;
    doa_theta   = linspace(-pi/2, pi/2, N_theta);       % general DOA vector 
    D_i         = D*linspace(1/N_D, 1, N_D);
    M_i         = 2:2+N_M-1;
    f_i         = f*linspace(1/N_F, 1, N_F);
    G_thetaM    = zeros(N_theta, N_M);
    G_thetaD    = zeros(N_theta, N_D);
    G_thetaF    = zeros(N_theta, N_F);
    
    for th = 1:N_theta
        for i = 1:N_D            
            m           = (0:M-1)';
            aULA_theta      = exp(2j*pi*D_i(i)/lambda*sin(doa_theta(th) - doa)*m); % general ULA vector       
            G_thetaD(th, i) = w'*aULA_theta;
        end
        
        for i = 1:N_M
            m           = (0:M_i(i)-1)';
            w_i = window(fhandle, M_i(i));
            aULA_theta      = exp(2j*pi*D/lambda*sin(doa_theta(th) - doa)*m); % general ULA vector       
            G_thetaM(th, i) = w_i'*aULA_theta;
        end
        
        for i = 1:N_F
            m           = (0:M-1)';
            aULA_theta      = exp(2j*pi*D/(c/f_i(i))*sin(doa_theta(th) - doa)*m); % general ULA vector       
            G_thetaF(th, i) = w'*aULA_theta;
        end
    end 
    
%% Display figures:
    
    figure('units','normalized','outerposition',[0 0 1 1])
    colordef white
    set(gcf, 'color', [1 1 1])
    fig = waterfall(doa_theta*180/pi, D_i*100, abs(G_thetaD'));
    set(fig, 'FaceColor', 'flat');
    set(fig, 'EdgeColor', 'k');
    set(fig, 'FaceVertexCData', [0.5 0.5 0.5])
    axis([-90 90 min(D_i)*100 max(D_i)*100])
    view(165, 65);
    xlabel 'DOA (degrees)'
    ylabel 'Separation (cm)'
    title(['Beam pattern vs. Sensor separation (M = 16, f = ' num2str(f) ' Hz)'])
    
    figure('units','normalized','outerposition',[0 0 1 1])
    colordef white
    set(gcf, 'color', [1 1 1])
    fig = waterfall(doa_theta*180/pi, M_i, abs(G_thetaM'));
    set(fig, 'FaceColor', 'flat');
    set(fig, 'EdgeColor', 'k');
    set(fig, 'FaceVertexCData', [0.5 0.5 0.5])
    axis([-90 90 min(M_i) max(M_i)])
    view(165, 65);
    xlabel 'DOA (degrees)'
    ylabel 'Number of sensors'
    title(['Beam pattern vs. Number of sensors (f = ' num2str(f) ' Hz)'])
    
    figure('units','normalized','outerposition',[0 0 1 1])
    colordef white
    set(gcf, 'color', [1 1 1])
    fig = waterfall(doa_theta*180/pi, f_i, abs(G_thetaF'));
    set(fig, 'FaceColor', 'flat');
    set(fig, 'EdgeColor', 'k');
    set(fig, 'FaceVertexCData', [0.5 0.5 0.5])
    axis([-90 90 min(f_i) max(f_i)])
    view(165, 65);
    xlabel 'DOA (degrees)'
    ylabel 'Frequency (Hz)'
    title(['Beam pattern vs. Frequency (M = 16, f = ' num2str(f) ' Hz)'])