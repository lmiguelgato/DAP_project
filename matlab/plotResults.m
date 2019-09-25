close all
clear
clc

addpath('../output')

max_num_sources = 1;

fileID  = fopen(['debug_' num2str(max_num_sources) 'data.txt']);

format = '%d:';
for i = 1:max_num_sources
    format = strcat(format, ' %f');
end

A       = textscan(fileID, format, 'Delimiter',',','EmptyValue',NaN);

meanDOAs = zeros(max_num_sources,1);
stdevDOAs = meanDOAs;

figure
for i = 1:max_num_sources
    meanDOAs(i)  = nanmean(A{i+1});
    stdevDOAs(i) = sqrt(nanvar(A{i+1}));
    
    hold on
    plot(A{1}, A{i+1},'LineWidth',2)
end
ylim([-184 184])
grid on

figure
for j = 1:max_num_sources
    polarplot(exp(1i*[meanDOAs(j) - stdevDOAs(j) meanDOAs(j) meanDOAs(j) + stdevDOAs(j)]/180*pi),':','LineWidth',4);
    hold on;
    polarplot(0.1*exp(1i*[-60 60 180 -60]/180*pi),'k','LineWidth',2);
    hold on;
end
    
for j = 1:max_num_sources
    hold on;
    polarplot([0.9*exp(1i*(meanDOAs(j) - stdevDOAs(j))/180*pi) 1.1*exp(1i*(meanDOAs(j) - stdevDOAs(j))/180*pi)],'k');
    hold on;
    polarplot([0.9*exp(1i*(meanDOAs(j) + stdevDOAs(j))/180*pi) 1.1*exp(1i*(meanDOAs(j) + stdevDOAs(j))/180*pi)],'k');
%     hold on;
%     polarplot(exp(1i*(meanDOAs(j) - stdevDOAs(j))/180*pi),'kx','LineWidth',5);
end

ax = gca;
ax.ThetaLim = [-180 180];
ax.RTickLabel = {''};
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';

