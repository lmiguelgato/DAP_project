close all
clear
clc

K = 20;

tau = 2048/44100;
acc = 1;

P_ = zeros(4);

x_ = [1;1;0;0];
y_all = zeros(2, K);
x_all = zeros(4, K);
for k = 1:K
    y = k*ones(2,1) + sqrt(40)*rand(2,1);
    [ x, P ] = kalmanFilter( x_, y, P_, acc, tau );
    y_all(:,k) = y;
    x_all(:,k) = x;
    x_ = x;
    P_ = P;
end

figure
plot(y_all(1,:), y_all(2,:))
hold on
plot(x_all(1,:), x_all(2,:))
grid on