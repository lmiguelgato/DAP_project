function [ x, P ] = kalmanFilter( x_, y, P, acc, tau )
% Kalman filter

% input:    x_  (previous state)
%           y   (current measurement)
%           P   (previous-state error posterior covariance matrix)
%           acc (acceleration from the movement model)
%           tau (observation period)

% output:   x   (current state)
%           P   (current-state error posterior covariance matrix)

% lmiguelgato@gmail.com

nvar = 0.001;   % measurement noise variance

% state transition model:
F = [1,   0, tau,   0;  % 4x4
     0,   1,   0, tau;
     0,   0,   1,   0;
     0,   0,   0,   1];

G = [tau^2/2,       0;  % 4x2
           0, tau^2/2;
         tau,       0;
            0,    tau];

Qw = sqrt([acc,   0;    % 2x2
             0, acc]);

Q = G*Qw*G';        % 4x4

% measurement model:
H = [1, 0, 0, 0;    % 2x4
     0, 1, 0, 0];

R = [nvar,    0;    % 2x2
        0, nvar];

%% prediction:
x_ = F*x_;          % 4x1
P_ = F*P*F' + Q;    % 4x4

%% update:
K = (P_*H')/(H*P_*H' + R);  % 4x2
x = x_ + K*(y - H*x_);      % 4x1
P = (eye(4) - K*H)*P_;      % 4x4

end