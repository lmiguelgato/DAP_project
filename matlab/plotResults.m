close all
clear
clc

addpath('../output')

max_num_sources = 2;

if max_num_sources == 1
    fileID  = fopen(['track_' num2str(max_num_sources) '_source.txt']);
else
    fileID  = fopen(['track_' num2str(max_num_sources) '_sources.txt']);
end

format = '%d:';
for i = 1:max_num_sources
    format = strcat(format, ' %f');
end

A = textscan(fileID, format, 'Delimiter',',','EmptyValue',NaN);


for i = 1:max_num_sources
    temp = A{1+i};
    for j = 1:length(temp)
        if temp(j) == 181
            temp(j) = NaN;
        end
    end
    A{1+i} = temp;
end

meanDOAs = zeros(max_num_sources,1);
stdevDOAs = meanDOAs;

figure
for i = 1:max_num_sources
    meanDOAs(i)  = nanmean(A{i+1});
    stdevDOAs(i) = sqrt(nanvar(A{i+1}));
    
    hold on
    plot(A{1}, A{i+1},'LineWidth',2)
end
ylim([-180.5 180.5])
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
end

ax = gca;
ax.ThetaLim = [-180 180];
ax.RTickLabel = {''};
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';
