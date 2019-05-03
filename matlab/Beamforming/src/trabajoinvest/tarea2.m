% Homework 2
% Only-sums Beamformer

    close all
    clear
    clc

    addpath('..\tools')
    
%% Parameters:

    fs      = 48e3;                         % sampling rate (Hz)
    f0      = 1.0e3;                        % 1st carrier frequency (narrow-band signal)
    f1      = 0.1e3;                        % 2nd carrier frequency (narrow-band signal)
    f       = [f0, f1];                     % vector of carrier frequencies

    N       = 1024;                         % signal duration
    t       = (0:N-1)'/fs;                  % time vector

    s0      = cos(2*pi*f(1)*t + rand*pi);   % 1st signal of interest
    s1      = cos(2*pi*f(2)*t + rand*pi);   % 2nd signal of interest
    
    k = 1;                          % index of signal to be separated

    c       = soundSpeed();         % speed of sound (m/s, t = 20°C by default)
    lambda  = c./f;                 % wavelength (meters)

    theta0  = 20;                   % 1st direction of arrival (degrees)
    theta1  = 0;                    % 2nd direction of arrival (degrees)
    theta   = [theta0, theta1];     % vector of directions of arrival (degrees)
    doa     = theta * pi/180;       % vector of directions of arrival (radians)

    M       = 16;                   % number of sensors
    D       = lambda(k)/2;          % separation of consecutive sensors (meters)
                                    % for a standard ULA, the inter-element spacing is half a wavelength
                                
    d       = D/lambda(k);          % this ratio should be less than 0.5 (grating lobes may appear otherwise)
    
%% Uniform linear array(ULA):

    % Sensor outputs:
    X = zeros(M,N);
    X(1,:) = s0+s1;
    for m = 1:M-1
        X(m+1,:) = delayFreq(s0, (m*D/c)*sin(doa(1)), fs) + ...
            delayFreq(s1, (m*D/c)*sin(doa(2)), fs);
    end
    
    y = sum(X)/16;
    % y = X'*ones(16,1)/16;
    
%% Power vs. angle:
    m       = (0:M-1)';                         % indices of sensors
    K       = 1000;
    angle   = linspace(-pi/2, pi/2, K);

    P_all = zeros(K, 1);
    for i = 1:K
        aULA_i  = exp(2j*pi*d*sin(angle(i))*m);       % ULA steering vector
        wULA_i  = aULA_i/sqrt(aULA_i'*aULA_i);
        Y_i = wULA_i'*X/sqrt(M)*2;
        P_all(i) = abs(Y_i)*abs(Y_i')/N;
    end
    
%% Display figures:

    figure
    plot(angle*180/pi, P_all/max(P_all), 'LineWidth', 2)
    axis([-90 90 0 1])
    grid on
    xlabel 'DOA (degrees)'
    title 'Normalized spectrum'    
    
    A = max(max(abs(y)), max(abs(X(1,:))))*1.2;

    figure
    plot(t*1e3, y, 'LineWidth', 2)
    hold on;
    plot(t*1e3, X(1,:), 'LineWidth', 2)
    axis([t(1)*1e3 t(end)*1e3 -A A])
    title 'Beamformer output'
    xlabel 'Time (ms)'
    