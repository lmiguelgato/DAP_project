close all
clear
clc

addpath('../output')

N = 1;

A = textread(['tabbed' num2str(N) 'data.txt'], '', 'delimiter', ',', ... 
                'emptyvalue', NaN);

% meanDOAs  = nanmean(A(:, 1:end-1));
% stdevDOAs = sqrt(nanvar(A(:, 1:end-1)));
            
meanDOAs  = nanmean(A(end-51:end, 1:end-1));
stdevDOAs = sqrt(nanvar(A(end-51:end, 1:end-1)));
            
disp(['Mean DOAs:          ' num2str(meanDOAs)])
disp(['Standard deviation: ' num2str(stdevDOAs)])
            
createfigure(A)
axis([1 size(A,1) -184 184])

N = size(A,1);

figure
for j = 1:size(A,2)-1
    polarplot(exp(1i*meanDOAs(j)/180*pi),'*','LineWidth',5);
    hold on;
end
    
for j = 1:size(A,2)-1
    hold on;
    polarplot(exp(1i*(meanDOAs(j) + stdevDOAs(j))/180*pi),'kx');
    hold on;
    polarplot(exp(1i*(meanDOAs(j) - stdevDOAs(j))/180*pi),'kx');
end

ax = gca;
ax.ThetaLim = [-180 180];
ax.RTickLabel = {''};
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';

