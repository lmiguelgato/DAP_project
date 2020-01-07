close all
clear
clc

addpath('../output')

enableML = 0;

max_num_sources = 4;

if max_num_sources == 1
    filename = strcat('track_', num2str(max_num_sources), '_source');
    fileID  = fopen(strcat(filename, '.txt'));
else
    filename = strcat('track_', num2str(max_num_sources), '_sources');
    fileID  = fopen(strcat(filename, '.txt'));
end

format = '%d:';
for i = 1:max_num_sources
    format = strcat(format, ' %f');
end

A = textscan(fileID, format, 'Delimiter',',','EmptyValue',NaN);

fileID  = fopen(strcat(filename, 'Kalman', '.txt'));
format = '%d %f %f';
K = textscan(fileID, format,'EmptyValue',NaN);

instantaneousSource = K{1};
instantaneousDOA = K{2};
kalmanDOA = K{3};

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

numKalmanDOAs   = zeros(max_num_sources,1);

temp = K{1};
for i = 1:max_num_sources
    for j = 1:length(temp)
        if temp(j) == i-1
            numKalmanDOAs(i) = numKalmanDOAs(i) + 1;
        end
    end
end

meanDOAs  = zeros(max_num_sources,1);
kalmanDOAs  = zeros(max_num_sources,1);
stdevDOAs = meanDOAs;
kalmanStdevDOAs = kalmanDOAs;
percentDOAs   = numDOAs/sum(numDOAs)*100;
percentKalmanDOAs   = numKalmanDOAs/sum(numKalmanDOAs)*100;

for i = 1:max_num_sources
    idx = find(instantaneousSource == i-1);
    meanDOAs(i)  = nanmean(A{i+1});
    kalmanDOAs(i) = nanmean(kalmanDOA(idx));
    stdevDOAs(i) = sqrt(nanvar(A{i+1}));
    kalmanStdevDOAs(i) = sqrt(nanvar(kalmanDOA(idx)));
end

firstDetection = sum(numDOAs);
for i = 1:max_num_sources
    tmp = find(~isnan(A{1+i}), 1);
    if tmp < firstDetection
        firstDetection = tmp;
    end
end

if (~enableML) 
    figure
    axis([1 sum(numDOAs)+firstDetection-1 -180.5 180.5])
    grid on

    labels = {};

    for i = 1:max_num_sources    
        idx = find(instantaneousSource == i-1);
        hold on; plot(instantaneousDOA(idx),'LineWidth',2)
        labels{end+1} = strcat('Instantaneous_', num2str(i));
        hold on; plot(kalmanDOA(idx),'LineWidth',2)
        labels{end+1} = strcat('Kalman_', num2str(i));
    end
    legend(labels)
end

%% ---------- kalman polar
if (~enableML) 
    figure
    for j = 1:max_num_sources
        if ~isnan(kalmanDOAs(j)) && ~isnan(stdevDOAs(j))
            tmp = exp(1i*[kalmanDOAs(j) - kalmanStdevDOAs(j) kalmanDOAs(j) kalmanDOAs(j) + kalmanStdevDOAs(j)]/180*pi);
            polarplot(tmp,':','LineWidth',ceil(8*percentKalmanDOAs(j)/100));
            hold on;
            text(imag(tmp(2)),real(tmp(2)),[num2str(round(percentKalmanDOAs(j))) ' %'])
        end
    end
    polarplot(0.1*exp(1i*[-60 60 180 -60]/180*pi),'k','LineWidth',2);

    for j = 1:max_num_sources
        hold on;
        polarplot([0.9*exp(1i*(kalmanDOAs(j) - kalmanStdevDOAs(j))/180*pi) 1.1*exp(1i*(kalmanDOAs(j) - kalmanStdevDOAs(j))/180*pi)],'k');
        hold on;
        polarplot([0.9*exp(1i*(kalmanDOAs(j) + kalmanStdevDOAs(j))/180*pi) 1.1*exp(1i*(kalmanDOAs(j) + kalmanStdevDOAs(j))/180*pi)],'k');
    end

    ax = gca;
    ax.ThetaLim = [-180 180];
    ax.RTickLabel = {''};
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
end

numDetectedSources = sum(percentKalmanDOAs > max(percentKalmanDOAs/2));

clc
if (~enableML) 
    disp(['Total number of sources detected: ' num2str(sum(percentKalmanDOAs>0)) ', found at: '])
    disp(num2str(kalmanDOAs(find(percentKalmanDOAs))))
    disp(' ')
end
disp(['Maximum-likelihood number of soures: ' num2str(numDetectedSources) ', found at: '])
disp(num2str(kalmanDOAs(find(percentKalmanDOAs > max(percentKalmanDOAs)/2))))


%% ---------- ML tracking

labels = {};
figure
MLindexes = find(percentKalmanDOAs > max(percentKalmanDOAs)/2);
for i = 1:length(MLindexes)
    idx = find(instantaneousSource == MLindexes(i)-1);
    hold on; plot(instantaneousDOA(idx),'LineWidth',2)
    labels{end+1} = strcat('Instantaneous_', num2str(i));
    hold on; plot(kalmanDOA(idx),'LineWidth',2)
    labels{end+1} = strcat('Kalman_', num2str(i));
end
legend(labels)

axis([1 sum(numDOAs)+firstDetection-1 -180.5 180.5])
grid on

%% ---------- ML kalman polar

figure
for idx = 1:length(MLindexes)
    j = MLindexes(idx);
    if ~isnan(kalmanDOAs(j)) && ~isnan(stdevDOAs(j))
        tmp = exp(1i*[kalmanDOAs(j) - kalmanStdevDOAs(j) kalmanDOAs(j) kalmanDOAs(j) + kalmanStdevDOAs(j)]/180*pi);
        polarplot(tmp,'LineWidth',ceil(8*percentKalmanDOAs(j)/100));
        hold on;
        text(imag(tmp(2)),real(tmp(2)),[num2str(round(percentKalmanDOAs(j))) ' %'])
    end
end
polarplot(0.1*exp(1i*[-60 60 180 -60]/180*pi),'k','LineWidth',2);
    
for idx = 1:length(MLindexes)
    j = MLindexes(idx);
    hold on;
    polarplot([0.9*exp(1i*(kalmanDOAs(j) - kalmanStdevDOAs(j))/180*pi) 1.1*exp(1i*(kalmanDOAs(j) - kalmanStdevDOAs(j))/180*pi)],'k');
    hold on;
    polarplot([0.9*exp(1i*(kalmanDOAs(j) + kalmanStdevDOAs(j))/180*pi) 1.1*exp(1i*(kalmanDOAs(j) + kalmanStdevDOAs(j))/180*pi)],'k');
end

ax = gca;
ax.ThetaLim = [-180 180];
ax.RTickLabel = {''};
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';