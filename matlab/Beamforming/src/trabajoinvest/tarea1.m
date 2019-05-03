% Homework 1
% Simulate radar detection

    close all
    clear
    clc

    addpath('..\tools')

%% Radar parameters:

    fs  = 48e3;         % sampling rate (Hz)
    dT  = 300e-3;       % repeat period (s)
    T   = 2.0;          % observation period (s)    

    % Radar signal (chirp)        
    tau = 20e-3;        % signal duration (s)
    N   = tau*fs;       % number of signal samples
    n   = (0:N-1)';     % discrete time    
    f1  = 5e3;          % initial frequency (Hz)
    f2  = 20e3;         % final frequency (Hz)

    % Instantaneous angle (radians)
    phi_inst = 2*pi*(f1*n + (f2-f1)*n.^2/2/N)/fs;  

    sn = cos(phi_inst);

%% Target parameters:

    L1  = 120;          % initial distance (m)
    L2  = 100;          % final distance (m)
    phi = 1 * (pi/180); % change in direction of arrival (radians)
    
    % Displacement (m)
    dL  = sqrt( L1^2 + L2^2 - 2*L1*L2*cos(phi) );

%% Simulate reflection:

    SNR = 0;            % signal-to-noise ratio (dB)
    c   = soundSpeed(); % sound speed (m/s, t = 20°C by default)
    
    t1 = 2*L1/c;        % 1st round trip time (s)
    t2 = 2*L2/c;        % 2nd round trip time (s)
    
    n1 = round(t1*fs);  n2 = round(t2*fs);

    % Reception
    nVar = var(sn/(L1^2))/10^(SNR/10);
    
    xn1 = sqrt(nVar)*randn(T*fs, 1);   xn1(n1:n1+N-1) = xn1(n1:n1+N-1) + sn/(L1^2);
    xn2 = sqrt(nVar)*randn(T*fs, 1);   xn2(n2:n2+N-1) = xn2(n2:n2+N-1) + sn/(L2^2);

%% Simulate reception:

    % Correlation receiver
    sn1 = xcorr(xn1, sn);
    sn2 = xcorr(xn2, sn);

    [m1, id1] = max(sn1);
    [m2, id2] = max(sn2);

    % Estimation of round trip times
    t_est1 = (id1 - T*fs + 1)/fs;
    t_est2 = (id2 - T*fs + 1)/fs;

    % Estimation of target distance
    L_est1 = t_est1*c/2;
    L_est2 = t_est2*c/2;

    % Estimation of displacement
    dL_est = sqrt( L_est1^2 + L_est2^2 - 2*L_est1*L_est2*cos(phi) );

    % Estimation of target speed
    v_est = dL_est/dT;  

    disp(['Estimated speed: ' num2str(v_est*3.6) ' km/h'])

%% Display figures:

    figure
    spectrogram(sn,64,60,128,fs,'yaxis')
    title(['Radar signal (chirp from ' num2str(f1/1e3) ...
        ' kHz to ' num2str(f2/1e3) ' kHz)'])

    
    t   = (1:T*fs)'/fs;         % sampled time (s)
    A = max(max(abs(xn1)), max(abs(xn2)))*1.2;
    
    figure, 
    subplot(211), plot(t, xn1)
    axis([t(1) t(end) -A A])
    title 'Received signals'
    ylabel '1st pulse'
    xlabel 'Time (s)'
    subplot(212), plot(t, xn2)
    axis([t(1) t(end) -A A])
    ylabel '2nd pulse'
    xlabel 'Time (s)'
    
    
    figure
    plot(t, sn1(T*fs:end))
    hold on
    plot(t, sn2(T*fs:end))
    title 'Correlation output'
    xlabel 'Time (s)'














