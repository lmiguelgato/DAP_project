close all
clear
clc

addpath('./tools')
addpath('./tools/bss_eval')

N = 3;

INpath = ['../corpus/corpus48000/clean-' num2str(N) 'source090180/'];
OUTpath = '../output/';

fs = 48e3;

T = 20;

L1 = 800000;
L2 = 950000;

D = 0.050*fs;
D1 = 0.010*fs;

L = round(T*fs);

[x, ~] = audioread([INpath 'wav_mic1.wav']);

S = zeros(N, L-D);
for i = 1:N
    [s, ~] = audioread([INpath 'pristine_channel' num2str(i) '.wav']);
    S(i, :) = s(1:L-D)';
end

for k1 = 1:N
    [y, ~] = audioread([OUTpath 'audio_' num2str(N) '_(' num2str(k1) ').wav']);
    y = y(D+1:L);
    
    for k2 = 1:N

    %% Evaluation of the separation algorithm:
    [s_target, e_interf, e_artif] = bss_decomp_filt(y(L1:L2)', k2, S(:, L1:L2), D1);

    [SDR, SIR, SAR] = bss_crit(s_target, e_interf, e_artif);

    [is_target, ie_interf, ie_artif] = bss_decomp_filt(x(L1:L2)', k2, S(:, L1:L2), D1);

    [SDRi, SIRi, SARi] = bss_crit(is_target, ie_interf, ie_artif);

    SIRd = SIR - SIRi;

    disp([num2str(k1) ', ' num2str(k2)])
    disp(['Input Signal to Interference Ratio:  ' num2str(SIRi) ' dB'])
    disp(['Output Signal to Interference Ratio: ' num2str(SIR) ' dB'])
    disp(['Improvement in Signal to Interference Ratio: ' num2str(SIRd) ' dB'])
    end
end

%     % Local evaluation:
%     [SDR_local, SIR_local, SAR_local] = ...
%         bss_crit(s_target, e_interf, zeros(1, L2-L1+1), e_artif, hanning(2048)', 1024);
%     % Local evaluation:
%     [iSDR_local, iSIR_local, iSAR_local] = ...
%         bss_crit(is_target, ie_interf, zeros(1, L2-L1+1), ie_artif, hanning(2048)', 1024);

% t       = (0:L2-L1)/fs;
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% plot(t, s_target, 'LineWidth', 2)
% hold on;
% plot(t, e_interf, 'LineWidth', 2)
% legend('s_{target}', 'e_{interf}')
% xlabel 'Time (s)'
% title 'Signal decomposition'
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% plot(t, is_target, 'LineWidth', 2)
% hold on;
% plot(t, ie_interf, 'LineWidth', 2)
% legend('is_{target}', 'ie_{interf}')
% xlabel 'Time (s)'
% title 'Signal decomposition'
% 
% figure('units','normalized','outerposition',[0 0 1 1])
% plot(SIR_local, 'LineWidth', 2)
% hold on
% plot(iSIR_local, 'LineWidth', 2)
% legend('SIR(t)')
% title 'Performance vs. time'

% player = audioplayer(y*40, fs);
% play(player)

% player = audioplayer(y, fs);
% play(player)
% 
% [x, ~] = audioread('wav_mic1.wav');
% 
% player = audioplayer(x, fs);
% play(player)

