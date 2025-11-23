function [] = plotPersist(frac, groups, groupNames, stimulusNames)
%%
figure
for g=1:length(groups)
    for k=1:length(stimulusNames)
        subplot(length(stimulusNames), length(groups), g+length(groups)*(k-1))
        data = frac(groups{g},:,k);
        hold on
        errorbar(mean(data), std(data) ./ sqrt(length(data)-1), 'r', 'Linewidth', 3)
        if k==1
            title(groupNames{g})
        elseif k==length(stimulusNames)
            xticks([1 2 3 4 5])
        end
        ylabel(stimulusNames{k})
        xlabel('Number of Sessions')
        ylim([0 .8])
    end
end
end