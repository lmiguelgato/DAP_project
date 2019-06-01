close all
clear
clc

N = 4;

addpath(['../output/clean-' num2str(N) 'source/'])

A = textread(['tabbed' num2str(N) 'data.txt'], '', 'delimiter', ',', ... 
                'emptyvalue', NaN);
				
if (size(A,1) > 50)
    meanDOAs  = nanmean(A(end-51:end, 1:end-1));
    stdevDOAs = sqrt(nanvar(A(end-51:end, 1:end-1)));
else
    meanDOAs  = nanmean(A(:, 1:end-1));
	stdevDOAs = sqrt(nanvar(A(:, 1:end-1)));
end
            
disp(['Mean DOAs:          ' num2str(meanDOAs)])
disp(['Standard deviation: ' num2str(stdevDOAs)])
            
createfigure(A)
axis([1 size(A,1) -184 184])

N = size(A,1);

figure
for j = 1:size(A,2)-1
    polarplot(exp(1i*meanDOAs(j)/180*pi),'*','LineWidth',10);
    hold on;
end
    
for j = 1:size(A,2)-1
    hold on;
    polarplot(exp(1i*(meanDOAs(j) + stdevDOAs(j))/180*pi),'kx','LineWidth',5);
    hold on;
    polarplot(exp(1i*(meanDOAs(j) - stdevDOAs(j))/180*pi),'kx','LineWidth',5);
end

ax = gca;
ax.ThetaLim = [-180 180];
ax.RTickLabel = {''};
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';

