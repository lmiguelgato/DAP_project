close all

A = textread('tabbed4data.txt', '', 'delimiter', ',', ... 
                'emptyvalue', NaN);

createfigure(A)
axis([1 size(A,1) -184 184])

N = size(A,1);

figure
for j = 1:size(A,2)-1
    %polarplot(exp(1i*A(N-10:N,j)/180*pi),'*','LineWidth',5);
    polarplot(exp(1i*nanmean(A(:,j))/180*pi),'*','LineWidth',5);
    ax = gca;
    ax.ThetaLim = [-180 180];
    ax.RTickLabel = {''};
    hold on;
end
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';

