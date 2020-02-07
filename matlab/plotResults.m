close all
clear
clc

folder = uigetdir('.', ...
         'Select directory where the location/tracking file are');

addpath(folder)

max_num_sources = 3;

if max_num_sources == 1
    filename = strcat('track_', num2str(max_num_sources), '_source');
    fileID  = fopen(strcat(filename, '.txt'));
else
    filename = strcat('track_', num2str(max_num_sources), '_sources');
    fileID  = fopen(strcat(filename, '.txt'));
end

if (fileID == -1)
    disp('Wrong number of sources. Aborting ...')
    return
end

format = '%d:';
for i = 1:max_num_sources
    format = strcat(format, ' %f');
end

A = textscan(fileID, format, 'Delimiter',',','EmptyValue',NaN);

numDOAs   = zeros(max_num_sources,1);

for i = 1:max_num_sources
    temp = A{1+i};
    for j = 1:length(temp)
        if temp(j) == 181
            temp(j) = NaN;
        else
            numDOAs(i) = numDOAs(i) + 1;
        end
    end
    A{1+i} = temp;
end

meanDOAs  = zeros(max_num_sources,1);
stdevDOAs = meanDOAs;
percentDOAs   = numDOAs/sum(numDOAs)*100;

figure
for i = 1:max_num_sources
    meanDOAs(i)  = nanmean(A{i+1});
    stdevDOAs(i) = sqrt(nanvar(A{i+1}));
    
    hold on
    plot(A{1}, A{i+1},'LineWidth',2)
end

firstDetection = sum(numDOAs);
for i = 1:max_num_sources
    tmp = find(~isnan(A{1+i}), 1);
    if tmp < firstDetection
        firstDetection = tmp;
    end
end

axis([1 sum(numDOAs)+firstDetection-1 -180.5 180.5])
grid on

fileID  = fopen(strcat(filename, 'Kalman', '.txt'));
format = '%d %f %f';
K = textscan(fileID, format,'EmptyValue',NaN);

instantaneousSource = K{1};
instantaneousDOA = K{2};
kalmanDOA = K{3};

labels = {};
for i = 1:max_num_sources    
    labels{end+1} = strcat('k-means_', num2str(i));
end

for i = 1:max_num_sources    
    idx = find(instantaneousSource == i-1);
    hold on; plot(instantaneousDOA(idx),'LineWidth',2)
    labels{end+1} = strcat('Instantaneous_', num2str(i));
    hold on; plot(kalmanDOA(idx),'LineWidth',2)
    labels{end+1} = strcat('Kalman_', num2str(i));
end
legend(labels)

figure
for j = 1:max_num_sources
    if ~isnan(meanDOAs(j)) && ~isnan(stdevDOAs(j))
        tmp = exp(1i*[meanDOAs(j) - stdevDOAs(j) meanDOAs(j) meanDOAs(j) + stdevDOAs(j)]/180*pi);
        polarplot(tmp,':','LineWidth',ceil(8*percentDOAs(j)/100));
        hold on;
        text(imag(tmp(2)),real(tmp(2)),[num2str(round(percentDOAs(j))) ' %'])
    end
end
polarplot(0.1*exp(1i*[-60 60 180 -60]/180*pi),'k','LineWidth',2);
    
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
