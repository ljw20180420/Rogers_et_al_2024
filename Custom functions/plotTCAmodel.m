function [] = plotTCAmodel(timeFactor,neuronFactor,trialFactor,colorCode)
[~,I] = sort(neuronFactor(:,1));
nDims = size(neuronFactor,2);
nTrials = size(trialFactor,1);
figure
for n=1:nDims
    neus = neuronFactor(:,n);

    subplot(nDims,3,1+(n-1)*3)
    plot(timeFactor(:,n),'Linewidth',3) %for relative time loadings, timeFactor(:,n)/sum(timeFactor(:,n))
    ylabel(n)
    ylim([0 max(timeFactor(:,n))+.2])
    xticks([0 10 35 55])
    ylim([0 max(timeFactor(:,n)+.2)])

    subplot(nDims,3,2+(n-1)*3)
    bar(neus(I),'c')

    subplot(nDims,3,3+(n-1)*3)
    scatter(1:nTrials, trialFactor(:,n), 50, colorCode, 'filled')
    hold on
    xline(8)
    hold on
    xline(16)
    hold on
    xline(22)
    hold on
    xline(28)
end
end