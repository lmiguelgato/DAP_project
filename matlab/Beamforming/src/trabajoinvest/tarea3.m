% Homework 3
% Weight-and-sum Beamformer

    close all
    clear
    clc

    addpath('..\tools')
    addpath('..\tools\bss_eval')
    
%% Parameters:

    winSel  = 8;

    fs      = 48e3;                         % sampling rate (Hz)
    f1      = 1.0e3;                        % 1st carrier frequency (narrow-band signal)
    f2      = 0.1e3;                        % 2nd carrier frequency (narrow-band signal)
    f       = [f1, f2];                     % vector of carrier frequencies

    N       = 2048;                         % signal duration
    t       = (0:N-1)/fs;                   % time vector

    s1      = cos(2*pi*f(1)*t + rand*pi);   % 1st signal of interest
    s2      = cos(2*pi*f(2)*t + rand*pi);   % 2nd signal of interest
    
    k       = 2;                            % index of signal to be separated

    c       = soundSpeed();         % speed of sound (m/s, t = 20°C by default)
    lambda  = c./f;                 % wavelength (meters)

    theta1  = 20;                   % 1st direction of arrival (degrees)
    theta2  = 0;                    % 2nd direction of arrival (degrees)
    theta   = [theta1, theta2];     % vector of directions of arrival (degrees)
    doa     = theta * pi/180;       % vector of directions of arrival (radians)

    M       = 16;                   % number of sensors
    D       = lambda(k)/2;          % separation of consecutive sensors (meters)
                                    % for a standard ULA, the inter-element spacing is half a wavelength
                                
    d       = D/lambda(k);          % this ratio should be less than 0.5 (grating lobes may appear otherwise)
    
%% Uniform linear array(ULA):

    % Sensor outputs:
    X = zeros(M,N);
    X(1,:) = s1+s2;
    for m = 1:M-1
        X(m+1,:) = delayTime(s1, (m*D/c)*sin(doa(1)), fs) + ...
            delayTime(s2, (m*D/c)*sin(doa(2)), fs);
    end
    
    switch winSel
        case 1
            fhandle = @rectwin;
        case 2
            fhandle = @bartlett;
        case 3
            fhandle = @hamming;
        case 4
            fhandle = @blackman;
        case 5
            fhandle = @blackmanharris;
        case 6
            fhandle = @chebwin;
        case 7
            fhandle = @gausswin;
        case 8
            fhandle = @taylorwin;
    end
            
    
    w = window(fhandle, M);
    
    y = w'*X/M;
    
%% Beam pattern:
    N_theta     = 1000;
    doa_theta   = linspace(-pi/2, pi/2, N_theta);       % general DOA vector 
    G_theta     = zeros(N_theta, 1);
    m           = (0:M-1)';                             % indices of sensors
    for th = 1:N_theta
        aULA_theta = exp(2j*pi*d*sin(doa_theta(th))*m); % general ULA vector
        G_theta(th) = w'*aULA_theta;
    end    
    
%% Evaluation of the separation algorithm:
    
    S = [s1; s2];

    [s_target, e_interf, e_artif] = bss_decomp_gain(y, k, S);

    [SDR, SIR, SAR] = bss_crit(s_target, e_interf, e_artif);

    [~, SIRi, ~] = bss_crit(S(k,:), sum(S,1)-S(k,:), zeros(1, N));

    SDRd = SDR - SIRi;
    SIRd = SIR - SIRi;

    disp(' ')
    disp(['Input Signal to Distortion Ratio:  ' num2str(SIRi) ' dB'])
    disp(['Output Signal to Distortion Ratio: ' num2str(SDR) ' dB'])
    disp(['Improvement in Signal to Distortion Ratio: ' num2str(SDRd) ' dB'])
    disp(['Improvement in Signal to Interference Ratio: ' num2str(SIRd) ' dB'])
    
    % Local evaluation:
    [SDR_local, SIR_local, SAR_local] = ...
        bss_crit(s_target, e_interf, zeros(1, N), e_artif, hanning(128)', 32);
    
%% Display figures:
    
    A = max(max(abs(y)), max(abs(X(1,:))))*1.2;

    figure('units','normalized','outerposition',[0 0 1 1])
    plot(t*1e3, X(1,:), 'LineWidth', 2)
    hold on;    
    plot(t*1e3, y, 'LineWidth', 2)
    axis([t(1)*1e3 t(end)*1e3 -A A])
    legend('s_{in}', 's_{out}')
    xlabel 'Time (ms)'
    title 'Beamformer input/output'
    
    figure('units','normalized','outerposition',[0 0 1 1])
    plot(t*1e3, s_target, 'LineWidth', 2)
    hold on;
    plot(t*1e3, e_interf, 'LineWidth', 2)
    hold on;
    plot(t*1e3, e_artif, 'LineWidth', 2)
    legend('s_{target}', 'e_{interf}', 'e_{artif}')
    A = 1.25*max([max(abs(X(1,:))), max(abs(s_target)), max(abs(e_interf))]);
    axis([min(t*1e3) max(t*1e3) -A A])
    xlabel 'Time (ms)'
    title 'Signal decomposition'
    
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
    
%     figure
%     plot(SDR_local, 'LineWidth', 2)
%     hold on;
%     plot(SIR_local, 'LineWidth', 2)
%     axis([1 length(SDR_local) 0 1.25*max(SIR_local)])
%     legend('SDR(t)', 'SIR(t)')
%     title 'Performance vs. time'
    
    
    
    
    
    
    