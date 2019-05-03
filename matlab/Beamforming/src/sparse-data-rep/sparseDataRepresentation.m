% Sparse data representation beamformer: 
% In this approach, the first step is to find a sparse representation of
% the array output data. Recently this approach has attracted many
% researchers' attention thanks to advances in the theory and methodology
% of sparse data reconstruction. It has a much better resolution than the
% conventional and MVDR beamformers at the expense of increased
% computational cost. Another advantage over the subspace methods is the
% improved robustness against signal coherence.

close all
clear
clc

addpath('..\tools')
addpath('..\tools\bss_eval')

%% Parameters:
fs      = 44.1e3;                       % sampling rate (Hz)
f0      = 1.0e3;                        % 1st carrier frequency (narrow-band signal)
f1      = 0.1e3;                        % 2nd carrier frequency (narrow-band signal)
f       = [f0, f1];                     % vector of carrier frequencies

N       = 1024;                         % signal duration
t       = (0:N-1)/fs;                   % time vector

s0      = cos(2*pi*f(1)*t + rand*pi);   % 1st signal of interest
s1      = cos(2*pi*f(2)*t + rand*pi);   % 2nd signal of interest
% s1 = trianglewave(3,N)*0.5;             % another type of signal

k = 1;                          % index of signal to be separated

c       = soundSpeed();         % speed of sound (m/s, t = 20°C by default)
lambda  = c./f;                 % wavelength (meters)

theta0  = 20;                   % 1st direction of arrival (degrees)
theta1  = 0;                    % 2nd direction of arrival (degrees)
theta   = [theta0, theta1];     % vector of directions of arrival (degrees)
doa     = theta * pi/180;       % vector of directions of arrival (radians)

M       = 12;                   % number of sensors
D       = lambda(k)/2;          % separation of consecutive sensors (meters)
                                % for a standard ULA, the inter-element spacing is half a wavelength
                                
d       = D/lambda(k);          % this ratio should be less than 0.5 (grating lobes may appear otherwise)
if d > 0.5
    disp(['Warning: d > 0.5, grating lobes may appear in the beampattern' ...
    'and create ambiguity during DOA estimation.'])
end