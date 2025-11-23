function [] = plot2metrics(AUCsMean, freezeFrac, sessionNames, groupNames, groupMembers, ylab1, ylab2, ylim1, ylim2)

figure
for n=1:length(groupNames)
    subplot(2, length(groupNames), n)
    barWithError(AUCsMean(:, groupMembers{n})');
    ylabel(ylab1)
    title(groupNames{n})
    ylim(ylim1)
    xticklabels(sessionNames)
end
for n=1:length(groupNames)
    subplot(2, length(groupNames), n + length(groupNames))
    barWithError(freezeFrac(:, groupMembers{n})');
    ylabel(ylab2)
    title(groupNames{n})
    ylim(ylim2)
    xticklabels(sessionNames)
end
end