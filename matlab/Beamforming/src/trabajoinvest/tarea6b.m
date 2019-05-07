% Homework 6
% Correlation-based DOA (two microphones, far-field assumption)
% Load sample audio (two channels)

    close all
    clear
    clc

    addpath('../tools')
    
%% Parameters:
    c       = soundSpeed(); % sound speed (m/s, t = 20�C by default)
    fc      = 1e3;          % emitter frequency
    d       = c/fc/2;       % separation between microphones
    L       = 1.5;          % emitter separation (m)
    
    N_theta = 20;       % number of theta samples
    theta_0 = -60;      % initial DOA
    theta_N = 60;       % final DOA
    theta   = linspace(theta_0/180*pi, theta_N/180*pi, N_theta);
    %theta   = linspace(theta_0/180*pi, theta_N/180*pi, N_theta) + 0.2*randn(1, N_theta);
    
    % Reception
    N_FFT   = 512;
    
    n1      = N_FFT/2;
    L1      = L;
    
    %% Load audio file
    [filename, path] = uigetfile( {'*.wav;','Audio files (*.wav)'; },'Pick channel 1');
    if isequal(filename,0)
       disp('The user cancelled the audio loading operation.')
       return;
    else
       [xn1,fs] = audioread(strcat((path),(filename)));      
    end
    
    %% Load audio file
    [filename, path] = uigetfile( {'*.wav;','Audio files (*.wav)'; },'Pick channel 2');
    if isequal(filename,0)
       disp('The user cancelled the audio loading operation.')
       return;
    else
       [xn2,fs] = audioread(strcat((path),(filename)));      
    end
    
    t       = (0:N_FFT*N_theta-1)'/fs;
    
    dt_est  = zeros(N_theta, 1);
    theta_est = dt_est;
    
    for step = 1:N_theta
        n1 = 2*N_FFT*(step-1)+1;
        n2 = 2*N_FFT*step;
        
        % Correlation receiver
        sn1 = xcorr(xn1(n1:n2), xn2(n1:n2));
        sn2 = ifft(   fft([xn1(n1:n2); zeros(2*N_FFT,1)]) ...
            .*conj(   fft([xn2(n1:n2); zeros(2*N_FFT,1)]) ));
        
        [m1, id1] = max(sn1);
        [m2, id2] = max(sn2);
        
        % Estimation of times
        t_est1 = (N_FFT*2 - id1)/fs;
        
        dt_est(step) = t_est1;        
        
        
        if abs(dt_est(step)) <= d/c
            theta_est(step) = asin(dt_est(step)*c/d);
        else
            disp(['Irrecoverable localization error at ' num2str(theta(step)*180/pi) '�'])
            if dt_est(step) < 0
                theta_est(step) = -pi/2;
            else
                theta_est(step) = pi/2;
            end
        end
        
        
    end

    figure1 = figure('units','normalized','outerposition',[0 0 1 1]);
    axes1 = axes('Parent',figure1);
    set(gcf, 'color', [1 1 1])
    plot(1:N_theta, (theta - theta_est')* 180/pi,'LineWidth',4)
    dTheta = max(abs(theta - theta_est'))* 180/pi;
    axis([1 N_theta -dTheta*1.5 dTheta*1.5])
    grid on
    xlabel 'Iteration'
    ylabel 'DOA error (degrees)'
    set(axes1,'FontSize',16);
    
    figure2 = figure('units','normalized','outerposition',[0 0 1 1]);
    axes2 = axes('Parent',figure2);
    set(gcf, 'color', [1 1 1])
    hold(axes2,'on');
    plot1 = plot([1:N_theta; 1:N_theta]',[theta', theta_est]* 180/pi,'LineWidth',4);
    set(plot1(2),'MarkerFaceColor',[0 0 0],'MarkerEdgeColor',[0 0 0],...
    'MarkerSize',8,'Marker','o','LineStyle',':','Color',[0 0 0]);
    axis([1 N_theta -90 90])
    box(axes2,'on');
    grid(axes2,'on');
    xlabel 'Iteration'
    ylabel 'DOA (degrees)'
    set(axes2,'FontSize',16);
    
    disp(' ')
    disp(['DOA error: ' num2str(rms((theta - theta_est'))* 180/pi) 'º (RMS)'])
    
