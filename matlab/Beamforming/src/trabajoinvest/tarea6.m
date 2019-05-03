% Homework 6
% Correlation-based DOA (two microphones, far-field assumption)

    close all
    clear
    clc

    addpath('..\tools')
    
%% Parameters:
    SNR     = 10;           % signal-to-noise ratio (dB)
    
    fs      = 48e3;         % sampling rate (Hz)
    fc      = 1e3;          % emitter frequency
    c       = soundSpeed(); % sound speed (m/s, t = 20°C by default)
    d       = c/fc/2;       % separation between microphones
    L       = 1.5;          % emitter separation (m)
    
    N_theta = 20;       % number of theta samples
    theta_0 = -60;      % initial DOA
    theta_N = 60;       % final DOA
    theta   = linspace(theta_0/180*pi, theta_N/180*pi, N_theta);
    %theta   = linspace(theta_0/180*pi, theta_N/180*pi, N_theta) + 0.2*randn(1, N_theta);
    
    % Reception
    N_FFT   = 1024;
    t       = (0:N_FFT*N_theta-1)'/fs;
    sn      = cos(2*pi*fc*t(1:N_FFT)) .* hamming(N_FFT);
    
    n1      = N_FFT/2;
    L1      = L;
    nVar    = var(sn/(L1^2))/10^(SNR/10);
    
    xn1     = sqrt(nVar)*randn(N_FFT*2, N_theta);         
    xn2     = sqrt(nVar)*randn(N_FFT*2, N_theta);   
    
    dt_est  = zeros(N_theta, 1);
    dt      = d*sin(theta)/c;
    theta_est = dt_est;
    
    for step = 1:N_theta
        n2 = n1 + round(dt(step)*fs);
        L2 = L1 - d*sin(theta(step));
        
        xn1(n1:n1+N_FFT-1, step) = xn1(n1:n1+N_FFT-1, step) + sn/(L1^2);
        xn2(n2:n2+N_FFT-1, step) = xn2(n2:n2+N_FFT-1, step) + sn/(L2^2);
        
        % Correlation receiver
        sn1 = xcorr(xn1(:, step), sn);
        sn2 = xcorr(xn2(:, step), sn);
        
        [m1, id1] = max(sn1);
        [m2, id2] = max(sn2);
        
        % Estimation of times
        t_est1 = (id1 - N_FFT + 1)/fs;
        t_est2 = (id2 - N_FFT + 1)/fs;
        
        dt_est(step) = t_est2 - t_est1;        
        
        
        if abs(dt_est(step)) <= d/c
            theta_est(step) = asin(dt_est(step)*c/d);
        else
            disp(['Irrecoverable localization error at ' num2str(theta(step)*180/pi) '°'])
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
    disp(['DOA error: ' num2str(rms((theta - theta_est'))* 180/pi) '° (RMS)'])
    