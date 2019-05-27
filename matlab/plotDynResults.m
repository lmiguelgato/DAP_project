close all
clear
clc

addpath('../output')

N = 2;

A = textread(['tabbed' num2str(N) 'data.txt'], '', 'delimiter', ',', ... 
                'emptyvalue', NaN);
            
createfigure(A)
axis([1 size(A,1) -184 184])

if N == 1
    [samples, dynDOA1] = textread('goldstandard.mdoa', ...
'Sample: %d\n\nDOAS:\n\n%f\n\n---\n\n');
hold on
plot(resample(dynDOA1, size(A,1), length(dynDOA1)), 'r', 'LineWidth',3)
end

if N == 2
    [samples, dynDOA1, dynDOA2] = textread('goldstandard.mdoa', ...
'Sample: %d\n\nDOAS:\n\n%f\n%f\n\n---\n\n');
hold on
plot(resample(dynDOA1, size(A,1), length(dynDOA1)), 'r', 'LineWidth',3)
hold on
plot(resample(dynDOA2, size(A,1), length(dynDOA2)), 'g', 'LineWidth',3)
end

if N == 3
    [samples, dynDOA1, dynDOA2, dynDOA3] = textread('goldstandard.mdoa', ...
'Sample: %d\n\nDOAS:\n\n%f\n%f\n%f\n\n---\n\n');
hold on
plot(resample(dynDOA1, size(A,1), length(dynDOA1)), 'r', 'LineWidth',3)
hold on
plot(resample(dynDOA2, size(A,1), length(dynDOA2)), 'g', 'LineWidth',3)
hold on
plot(resample(dynDOA3, size(A,1), length(dynDOA3)), 'c', 'LineWidth',3)
end

if N == 4
    [samples, dynDOA1, dynDOA2, dynDOA3, dynDOA4] = textread('goldstandard.mdoa', ...
'Sample: %d\n\nDOAS:\n\n%f\n%f\n%f\n%f\n\n---\n\n');
hold on
plot(resample(dynDOA1, size(A,1), length(dynDOA1)), 'r', 'LineWidth',3)
hold on
plot(resample(dynDOA2, size(A,1), length(dynDOA2)), 'g', 'LineWidth',3)
hold on
plot(resample(dynDOA3, size(A,1), length(dynDOA3)), 'c', 'LineWidth',3)
hold on
plot(resample(dynDOA4, size(A,1), length(dynDOA4)), 'y', 'LineWidth',3)
end
