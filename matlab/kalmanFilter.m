function [ x, P ] = kalmanFilter( x_, y, P, acc, tau )
% Kalman filter
% input: x_ (previous state)
% output: x (current state)
% lmiguelgato@gmail.com

nvar = 0.001;

F = [1,   0, tau,   0; 
     0,   1,   0, tau;
     0,   0,   1,   0;
     0,   0,   0,   1];
 
G = [tau^2/2,       0;
           0, tau^2/2;
         tau,       0;
            0,    tau];
        
Qw = sqrt([acc,   0;
             0, acc]);

R = [nvar,    0;
        0, nvar];
    
Q = G*Qw*G';
 
H = [1, 0, 0, 0;
     0, 1, 0, 0];
 
%% prediction:
x_ = F*x_;          % 4x1
P_ = F*P*F' + Q;    % 4x4

%% update:
K = (P_*H')/(H*P_*H' + R);  % 4x2
x = x_ + K*(y - H*x_);      % 4x1
P = (eye(4) - K*H)*P_;

end

