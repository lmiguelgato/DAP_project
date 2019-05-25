close all
clear
clc

N = 2;

A = textread(['tabbed' num2str(N) 'data.txt'], '', 'delimiter', ',', ... 
                'emptyvalue', NaN);

if N == 1
    [samples, dynDOA1] = textread('goldstandard.mdoa', ...
'Sample: %d\n\nDOAS:\n\n%f\n\n---\n\n');
end

if N == 2
    [samples, dynDOA1, dynDOA2] = textread('goldstandard.mdoa', ...
'Sample: %d\n\nDOAS:\n\n%f\n%f\n\n---\n\n');
end

createfigure(A)
axis([1 size(A,1) -184 184])
hold on
plot(resample(dynDOA1, size(A,1), length(dynDOA1)), 'r', 'LineWidth',3)
hold on
plot(resample(dynDOA2, size(A,1), length(dynDOA1)), 'r', 'LineWidth',3)