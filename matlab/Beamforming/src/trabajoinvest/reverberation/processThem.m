close all
clear
clc

[y_exp1, fs] = audioread('EXP_local_abierto_personas_dentro.wav');
[y_exp2, ~] = audioread('EXP_local_abierto_personas_fuera.wav');
[y_exp3, ~] = audioread('EXP_local_cerrado_personas_dentro.wav');

[y_lin1, ~] = audioread('LIN_local_abierto_personas_dentro.wav');
[y_lin2, ~] = audioread('LIN_local_abierto_personas_fuera.wav');
[y_lin3, ~] = audioread('LIN_local_cerrado_personas_dentro.wav');

[y_mls1, ~] = audioread('MLS_local_abierto_personas_dentro.wav');
[y_mls2, ~] = audioread('MLS_local_abierto_personas_fuera.wav');
[y_mls3, ~] = audioread('MLS_local_cerrado_personas_dentro.wav');

N = length(y_exp1);
t = (0:N-1)/fs;

n_exp1 = y_exp1(N/2:end);
n_exp2 = y_exp2(N/2:end);
n_exp3 = y_exp3(N/2:end);

n_lin1 = y_lin1(N/2:end);
n_lin2 = y_lin2(N/2:end);
n_lin3 = y_lin3(N/2:end);

n_mls1 = y_mls1(N/2:end);
n_mls2 = y_mls2(N/2:end);
n_mls3 = y_mls3(N/2:end);

var_exp1 = var(n_exp1);
var_exp2 = var(n_exp2);
var_exp3 = var(n_exp3);

var_lin1 = var(n_lin1);
var_lin2 = var(n_lin2);
var_lin3 = var(n_lin3);

var_mls1 = var(n_mls1);
var_mls2 = var(n_mls2);
var_mls3 = var(n_mls3);

var_exp = [var_exp1, var_exp2, var_exp3]
var_lin = [var_lin1, var_lin2, var_lin3]
var_mls = [var_mls1, var_mls2, var_mls3]

figure
plot(t, y_exp1)



