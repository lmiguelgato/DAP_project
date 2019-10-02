close all
clear
clc

addpath('../output')

max_num_sources = 1;

INpath = ['../corpus/' num2str(max_num_sources) ' Source/'];

%% loading ground truth:
[mdoa.name, mdoa.path] = uigetfile('*.mdoa', ...
    'Select ground truth.', INpath);

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

createfigure(A)
ylim([-180.5 180.5])

mdoa.format = 'Sample: %d\n\nDOAS:\n\n';
for i = 1:max_num_sources
    mdoa.format = strcat(mdoa.format, '%f\n');
end
mdoa.format = strcat(mdoa.format, '\n---\n\n');

if max_num_sources == 1
    [samples, dynDOA1] = textread([mdoa.path, mdoa.name], mdoa.format);
    hold on
    plot(resample(dynDOA1, size(A,1), length(dynDOA1)), 'r', 'LineWidth',3)
end

if max_num_sources == 2
    [samples, dynDOA1, dynDOA2] = textread([mdoa.path, mdoa.name], mdoa.format);
    hold on
    plot(resample(dynDOA1, size(A,1), length(dynDOA1)), 'r', 'LineWidth',3)
    hold on
    plot(resample(dynDOA2, size(A,1), length(dynDOA2)), 'g', 'LineWidth',3)
end