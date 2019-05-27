close all
clear
clc

addpath('./tools')
addpath('./tools/bss_eval')
addpath('./corpus/corpus48000/clean-4source')

k = 1;

[s1, fs] = audioread('pristine_channel1.wav');
[s2, ~] = audioread('pristine_channel2.wav');
[s3, ~] = audioread('pristine_channel3.wav');
[s4, ~] = audioread('pristine_channel4.wav');

[y, ~] = audioread('audio_4_(2).wav');

N = min([length(s1), length(s2), length(s3), length(s4), length(y)]);
N = floor(N/1024)*1024;

%% Evaluation of the separation algorithm:
S = [s1(1:N)'; s2(1:N)'; s3(1:N)'; s4(1:N)'];

[s_target, e_interf, e_artif] = bss_decomp_gain(y(1:N)', k, S);

[SDR, SIR, SAR] = bss_crit(s_target, e_interf, e_artif);

[~, SIRi, ~] = bss_crit(S(k,:), sum(S,1)-S(k,:), zeros(1, N));

SIRd = SIR - SIRi;

disp(' ')
disp(['Input Signal to Interference Ratio:  ' num2str(SIRi) ' dB'])
disp(['Output Signal to Interference Ratio: ' num2str(SIR) ' dB'])
disp(['Improvement in Signal to Interference Ratio: ' num2str(SIRd) ' dB'])

% Local evaluation:
[SDR_local, SIR_local, SAR_local] = ...
    bss_crit(s_target, e_interf, zeros(1, N), e_artif, hanning(2048)', 1024);

t       = (0:N-1)/fs;

figure('units','normalized','outerposition',[0 0 1 1])
plot(t, s_target, 'LineWidth', 2)
hold on;
plot(t, e_interf, 'LineWidth', 2)
legend('s_{target}', 'e_{interf}')
xlabel 'Time (s)'
title 'Signal decomposition'

figure('units','normalized','outerposition',[0 0 1 1])
plot(SIR_local, 'LineWidth', 2)
legend('SIR(t)')
title 'Performance vs. time'

% player = audioplayer(s_target*40, fs);
% play(player)

% player = audioplayer(y, fs);
% play(player)
% 
% [x, ~] = audioread('wav_mic1.wav');
% 
% player = audioplayer(x, fs);
% play(player)

