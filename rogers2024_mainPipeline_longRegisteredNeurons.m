

%% In Vivo Imaging Analysis Pipeline for Rogers et al., 2024
%
%
%Each part depends on previous parts and will not run without one another
%unless specified.
%
%In the following sections you will:
% 1. Load the data of all animals in a given experimental group from
% storeRogers2024Data.m
% 2. Extract the traces of longitudinally registered cells from all files
% and align to user-defined stimulus times or times of interest
% 3. Plot spatial coordinates of longitudinally registered cells
% 4. Extract and plot average traces in terms of % baseline of all neurons 
% averaged over all trials in each session
% 5. Identify cells that are upregulated or downregulated in response to
% stimuli in a given session according to their average trace
% 6. Measure the number of sessions cells tend to maintain stimulus
% responsiveness
% 7. Calculate freezing encoding in each neuron in each session and compare
% within animals across groups
% 8. Align freezing to trials
% 9. Load & plot TCA models, identifying session-dominant components and
% correlating to freezing
% 10. Calculate and plot normalized trial factors.
% 11. Calculate and plot normalized temporal factors.
% 12. Assess cumulative distribution of neuron weights
% 13. Identify Acq-Dominant, Ext1-Dominant, and Ext3-Dominant ensembles and
% their overlaps
% 14. Use a Fisher linear decoder to predict group based on average
% 15. Calculate changed activity (z-score relative to Acquisition) of each ensemble over time 
% 16. Reconstruct example neurons from TCA
% 17. Stimulus response properties of TCA ensemble neurons.
%
% To run non-shock mice, change for loops that are a=1:25 to a=26:28
%
% For aesthetic purposes, most data was exported from these sections in csv
% and txt files and imported in Prism, where most statistical tests were
% also done
%
%Written by Sophie A. Rogers, Corder Laboratory, University of Pennsylvania




%% 1. Load the data of all animals in a given experimental group this will take a minute or so


%the following structure contains all the raw data necessary to run the following code.
%the structure is broken down by animals (structures).
%within each animal field are subfields for sessions (structure), 
%cell registration table across sessions (table), 
%cell spatial properties (table), and tca output(table)

%sessions has 5 subfields (6 if broken videos), one for each session
%(habituation, acquisition, extinction 1, extinction 2, extinction 3). in
%each session subfield is calcium, a table of traces; stimulusTimes, a
%vector of stimulus delivery times; and freezing, a table of freezing

addpath('D:\z_serious_works\Rogers_et_al_2024/Custom functions/')
addpath('D:\z_serious_works\Rogers_et_al_2024/Rogers2024_CalciumData_AllAnimals/mm39/')
load('D:\z_serious_works\Rogers_et_al_2024/rogers2025CalciumData.mat')

% rogers2024.trialStructure.trialStart. defined here as how much time back from stimulus times to begin collecting data. (seconds)
% rogers2024.trialStructure.trialEnd. how much time after stimulus start to stop collecting data. (seconds)
downSamplingRate = 15;
%length of trial (seconds)
trialLength = rogers2024.trialStructure.trialStart + rogers2024.trialStructure.trialEnd + 1;

%% 2. Extract traces of longitudinal cells and align to stimuli

%CUSTOM FUNCTIONS REQUIRED:
%   crossSeshAuto
%   aligntraces

for a=1:length(rogers2024.animals)
        
    %longitudinally register and spatially locate cells across 5 sessions
    [longRegistered, coordinates{a}] = crossSeshAuto( ...
        rogers2024.animals(a).cellreg, ...
        {rogers2024.animals(a).sessions.calcium}, ...
        rogers2024.animals(a).cellprops ...
    );
   
    %repair broken recordings in 3 animals
    if a == 12
        longRegistered = {
            longRegistered{1}, ...
            longRegistered{2}, ...
            [
                longRegistered{3};
                longRegistered{4}
            ], ...
            longRegistered{5}, ...
            longRegistered{6}
        };
    elseif a == 13
        longRegistered = {
            [
                longRegistered{1}; 
                longRegistered{1}(end-16.8*rogers2024.trialStructure.samplingRate/2+1:end,:);
                longRegistered{2}(1:16.8*rogers2024.trialStructure.samplingRate/2,:);
                longRegistered{2}
            ], ...
            longRegistered{3}, ...
            longRegistered{4}, ...
            longRegistered{5}, ...
            longRegistered{6}
        };
    elseif a == 28
        longRegistered = {
            [
                longRegistered{1};
                longRegistered{2}
            ], ...
            longRegistered{3}, ...
            longRegistered{4}, ...
            longRegistered{5}, ...
            longRegistered{6}
        };
    end
    
    %loop through sessions, downsample registered traces to 15Hz and store.
    %also store stimulus delivery times
    for m=1:rogers2024.trialStructure.nSessions
        lr{a,m} = resample(longRegistered{m}, downSamplingRate, rogers2024.trialStructure.samplingRate);
    end
    
    %store number of cells
    nCells(a) = size(longRegistered{1}, 2);
    
    %align traces to stimulus times and store 1Hz data at the predefined
    %trial windows
    poolMat{a} = aligntraces( ...
        longRegistered, ...
        {rogers2024.animals(a).sessions.stimulusTimes}, ...
        rogers2024.trialStructure.nSessions, ...
        rogers2024.trialStructure.samplingRate, ...
        rogers2024.trialStructure.trialStart, ...
        rogers2024.trialStructure.trialEnd ...
    );
    
    %transform into s x c x t tensor and store
    tensors{a} = zeros(trialLength, nCells(a), rogers2024.trialStructure.nTrialsTotal);
    for t=1:rogers2024.trialStructure.nTrialsTotal
        tensors{a}(:,:,t) = poolMat{a}(trialLength*(t-1)+1:trialLength*t,:);
    end
end
poolMat = [poolMat{:}];


%% 3.Plot spatial coordinates of longitudinally registered cells. Fig 2B 

%enter animal number of interest; dataset 1 animal 1 is the representative image
propM39 = {
    readtable('mm39hab-props.csv'), ...
    readtable('mm39acq-props.csv'), ...
    readtable('mm39ext1-props.csv'), ...
    readtable('mm39ext2-props.csv'), ...
    readtable('mm39ext3-props.csv')
};

%plot cells at x coordinate by y coordinate at their size in pixels, scaled up to be visible
colors = {[.5 .5 .5], [1 0 0], [0.9290 0.6940 0.1250], [0 1 0], [0 1 1]};
figure
for m=1:rogers2024.trialStructure.nSessions
    hold on
    scatter(propM39{m}.CentroidX, propM39{m}.CentroidY, propM39{m}.Size * (20-3*(m-1)), colors{m}, 'filled')
end
hold on
scatter(coordinates{1}(:,1), coordinates{1}(:,2), coordinates{1}(:,3) * 8, [0.4940 0.1840 0.5560]) 
%% 4. for fig 2D - heatmaps of average traces

%CUSTOM FUNCTIONS
% percBaselineAvg

%storage cell average session activity of each animal by session
%loop through animals
for a=1:length(rogers2024.animals) %can look at specific groups by replacing "1:length(rogers2024.animals)" with saline, responders, nonresponders, or non-shock
    %call t x c x T data from tensors and convert to time x cells matrix
    data = tensors{a};
    %zscore to baseline
    data = (data - mean(data(rogers2024.stimuli(1).times, :, :))) ./ std(data(rogers2024.stimuli(1).times, :, :));
    
    %loop through sessions
    for m=1:rogers2024.trialStructure.nSessions
        %store median zscore in time x cells x session tensor
        avPBt{a}{m} = median(data(:, :, rogers2024.sessions(m).trials), 3);
    end
    avPBt{a} = cat(3, avPBt{a}{:});
end

%pool cells from all animals within sessions
%loop through animals to pool mean trial activity in a session
for g=1:length(rogers2024.groups)
    for a = rogers2024.groups(g).members
        %loop through sessions
        for m=1:rogers2024.trialStructure.nSessions
            %store all average traces in a session from every animal
            poolAvgAct{m,g}{a} = avPBt{a}(:,:,m);
        end
    end
    for m=1:rogers2024.trialStructure.nSessions
        poolAvgAct{m,g} = [poolAvgAct{m,g}{:}];
    end
end

%plot
figure
for g=1:length(rogers2024.groups)
    %sort data to shock
    [~, I] = sort(mean(poolAvgAct{2,g}(rogers2024.stimuli(4).times, :)));
    length(I)
    for m=1:rogers2024.trialStructure.nSessions
        subplot(length(rogers2024.groups), rogers2024.trialStructure.nSessions, m + 5 * (g - 1))

        data2plot = poolAvgAct{m,g};
        h=heatmap(data2plot(:,I)','Colorlimits',[-2 2]);
        colormap(hot)
        if g==1
            title(rogers2024.sessions(m).name)
        end
        ylabel('Cells')
        if g==length(rogers2024.groups)
            xlabel('Time')
        end
        h.XDisplayLabels = repmat({''}, 1, size(h.ColorData, 2));
        h.YDisplayLabels = repmat({''}, 1, size(h.ColorData, 1));
        if m < rogers2024.trialStructure.nSessions
            h.ColorbarVisible = 'off';
        end
        grid off
    end
end



%% 5. Identify cells that are upregulated or downregulated in response to stimuli in a given session according to their average trace for fig 2J-M

%CUSTOM FUNCTIONS
%permTest shock = 56:60 instead of rogers2024.stimuli(4).times
stimuli = {
    rogers2024.stimuli(2).times, ...
    rogers2024.stimuli(3).times, ...
    [rogers2024.stimuli(2).times rogers2024.stimuli(3).times], ...
    56:60
};

for a=1:length(rogers2024.animals)
    for m=1:rogers2024.trialStructure.nSessions        
        %for tone, trace, tone+trace, shock
        for k=1:length(stimuli)
            sigCells{a,m,k} = isSig(avPBt{a}(:,:,m), rogers2024.stimuli(1).times, stimuli{k});
        end
    end
end

%calculate fractions of activated neurons
for a=1:length(rogers2024.animals)
    for m=1:rogers2024.trialStructure.nSessions
        for k=1:length(stimuli)
            fracs(a,m,k) = length(sigCells{a,m,k}{1}) / nCells(a);
        end
    end
end

%plot
figure
for k=1:length(stimuli)
    subplot(2, 2, k)
    for g = 1:length(rogers2024.groups)
        errorbar(mean(fracs(rogers2024.groups(g).members, :, k)), std(fracs(rogers2024.groups(g).members, :, k)) ./ sqrt(length(rogers2024.groups(g).members)))
        hold on
    end
end



%% 6. Measure the number of sessions cells tend to maintain stimulus responsiveness. For fig 2N-O

%CUSTOM FUNCTIONS
%permTest

%loop through animals to identify tone, trace, and shock responsive neurons
for a = 1:length(rogers2024.animals)
    %call average activity tensor 
    %loop through stimuli
    for k = 1:length(stimuli)
        counter = zeros(rogers2024.trialStructure.nSessions, nCells(a));
        for c=1:nCells(a)
            for m=1:rogers2024.trialStructure.nSessions
                counter(m,c) = ismember(c, sigCells{a,m,k}{1});           
            end
        end
        
        for m=1:rogers2024.trialStructure.nSessions
            frac(a,m,k) = sum(sum(counter)==m) ./ nCells(a);
        end
    end
    disp(a)
 end
    

%plot
plotPersist(frac, {rogers2024.groups.members}, {rogers2024.groups.name}, {'tone', 'trace', 'tone+trace', 'shock'})

%% 7. Measure the freezing encoding of single neurons. For fig 2P-Q, Supp Fig. 2B
%Representative image is animal g=10, session=3

%initialize storage matrix for the fraction of neurons and average freezing
%encoding

%loop through animals
for a = 1:length(rogers2024.animals)
    %loop through sessions
    for m = 1:rogers2024.trialStructure.nSessions
        data = lr{a,m};    % downsampled traces
        %remove temporal offset between miniscope and freezing recording (seconds), z-score over the session
        ac = zscore(data(rogers2024.trialStructure.recordingOffset * downSamplingRate:end, :));
        
        freezing = table2array(rogers2024.animals(a).sessions(m).freezing);      %loop through freezing
        freezing(freezing==100) = 1;  %turn into binary
        
        if length(freezing) > size(ac, 1)  %make matrices the same length
            freezing = freezing(1:size(ac, 1));
        else
            ac = ac(1:length(freezing), :);
        end
        
        %initialize storage vectors for auROC and cell indices for auROC > .6
        AUCs = [];
        freezeEncoding{a,m} = [];
        
        %loop through neurons
        for c=1:size(ac, 2)                
            %create shuffle vector
            randN = randperm(round(length(freezing)));
            
            %pick random 50% of activity freezing values for training
            %set
            trainAct = ac(randN(1:round(end / 2)), c);
            trainFreeze = freezing(randN(1:round(end / 2)));
            
            %pick random 50% of activity freezing values for test
            %set
            testAct = ac(randN(round(end / 2) + 1:end), c);
            testFreeze = freezing(randN(round(end / 2) + 1:end));
            
            %train logistic regressor
            [B, ~, STATS] = glmfit(trainAct, trainFreeze, 'binomial', 'Link', 'logit');
            
            %test
            pred = glmval(B, testAct, 'logit');
            
            %calculate auROC
            [~, ~, ~, AUCs(c)] = perfcurve(testFreeze, pred, 1, 'XCrit', 'FPR', 'YCrit', 'TPR');
           
            if STATS.p(2)<0.0001   
                %record their index
                freezeEncoding{a,m} = [freezeEncoding{a,m}, c];
            end
        end
        
        %record average auROC of animal in session
        AUCsMean(m,a) = mean(AUCs);
        
        % ac(isnan(ac)) = 0;        
        % 
        % [~, I] = sort(AUCs);
        % cofI = I(ismember(I, freezeEncoding{a,m}))';
        % 
        % %Uncomment for representative image
        % figure
        % subplot(211)
        % heatmap(freezing(1241:1440)')
        % grid off
        % subplot(212)
        % heatmap(zscore(ac(1241:1440, cofI))', 'colorlimits', [-2 2])
        % colormap(hot)
        % grid off
    end
  disp(a)
end

%plot
for a = 1:length(rogers2024.animals)
    blPop = freezeEncoding{a,1};
    for m=2:rogers2024.trialStructure.nSessions
        pop = freezeEncoding{a,m};
        pop(ismember(pop, blPop)) = [];
        freezeFrac(m-1,a) = length(pop) / nCells(a);
    end
end

plot2metrics( ...
    AUCsMean, ...
    freezeFrac, ...
    {rogers2024.sessions.name}, ...
    {rogers2024.groups.name}, ...
    {rogers2024.groups.members}, ...
    'Mean freezing encoding', ...
    'Fraction of Cells', ...
    [.5 .7], ...
    [0 .7] ...
)
%% 8. Extract freezing and align to stimuli
%calculate tensors from Williams et al. TCA code (https://github.com/ahwillia/tensortools)
% load tensors 

%get trial by trial freezing
%get standard stimulus times
for m=1:rogers2024.trialStructure.nSessions
    timesF{m} = rogers2024.animals(24).sessions(m).stimulusTimes - rogers2024.trialStructure.recordingOffset;
end

%create second by second % freezing
for a=1:length(rogers2024.animals)
    for m=1:rogers2024.trialStructure.nSessions
        %read each animal's binary freezing data from each session
        fr{m} = table2array(rogers2024.animals(a).sessions(m).freezing);
        
        %downsample to 1 second by thresholding (if freezing is greater
        %than 50% in that second
        fzH = [];
        for ii=1:round(length(fr{m}) / downSamplingRate) - 1
            val = mean(fr{m}((ii - 1) * downSamplingRate + 1:ii * downSamplingRate));
            if val>=50
                fzH(ii) = 1;
            else
                fzH(ii) = 0;
            end
        end
        
        %get that animal's stimulus times
        for t=1:length(rogers2024.animals(a).sessions(m).stimulusTimes)
            % 这好像不对
            stimTime = round( ...
                rogers2024.animals(a).sessions(m).stimulusTimes(t) - rogers2024.animals(a).sessions(m).stimulusTimes(1) + rogers2024.trialStructure.recordingOffset ...
            );
            fzHz{a,m}{t} = fzH(stimTime - rogers2024.trialStructure.trialStart:stimTime + rogers2024.trialStructure.trialEnd);
        end
        fzHz{a,m} = [fzHz{a,m}{:}];
    end

    %fix broken video
    timesf = timesF;
    if a==5
        timesf{3} = [120; 120+round(2844/15); 234+round(2844/15); 331+round(2844/15); 441+round(2844/15); 543+round(2844/15)];
    end
    
    %align freezing to stimulus times
    freezingInTrials{a} = aligntraces( ...
        fr, ...
        timesf, ...
        rogers2024.trialStructure.nSessions, ...
        downSamplingRate, ...
        rogers2024.trialStructure.trialStart, ...
        rogers2024.trialStructure.trialEnd ...
    );
end


%take freezing 10 seconds before each tone to 50 seconds after on each
%trial and concatenate
for a=1:length(rogers2024.animals)
    for t = 1:rogers2024.trialStructure.nTrialsTotal
       freTot(t, a) = mean(freezingInTrials{a}(trialLength * (t - 1) + 1:trialLength * t));
    end
end

%% 9. Load & plot TCA models, identifying session-dominant components and correlating to freezing for fig. 3A,D-F

nDims = 5;

%to plot representative TCAs, uncomment the figure file in the for-loop -
%the representative image in the paper is from animal 12

for a=1:length(rogers2024.animals)
    if a < 22 || a > 25
        tca = table2array(rogers2024.animals(a).tca(2:end, :));
    else
        tca = table2array(rogers2024.animals(a).tca(:, 2:end));
    end
    %store neuron factors - used in Fig. 3G, 4-5
    timFac{a} = tca(1:trialLength,:);
    neuFac{a} = tca(trialLength + 1:end - rogers2024.trialStructure.nTrialsTotal, :);
    triFac{a} = tca(end - rogers2024.trialStructure.nTrialsTotal + 1:end, :);
    
    %Identify session-dominant factors - used in Fig. 3D-F, 4, 5
    for m=1:rogers2024.trialStructure.nSessions
        facLoads = mean(triFac{a}(rogers2024.sessions(m).trials, :));
        [~, domFacs(a,m)] = max(facLoads);
    end
    
    %calculate strength of each factor in each session - used in Fig. D,E
    for m=1:rogers2024.trialStructure.nSessions
        facLoads = mean(triFac{a}(rogers2024.sessions(m).trials, :));
        % 这好像不对
        strength(a,m) = facLoads(domFacs(a,m))/sum(facLoads(domFacs(a,:)));
    end

    %uncomment to plot representative figure
    %plotTCAmodel(timFac{a}, neuFac{a}, triFac{a}, freTot(:,a))
end

for m=1:rogers2024.trialStructure.nSessions
    for n=1:rogers2024.trialStructure.nSessions
        [r(n,m), p(n,m)] = corr(mean(freTot(rogers2024.sessions(n).trials, :))', strength(:, m));
    end
end
r(p>0.05) = 0;
figure
heatmap(r)
colormap(hot)

%% 10. Plot normalized trial factors. for Supp Fig. 4A-E

for a=1:length(rogers2024.animals)
    %call trial weights from a given animal    
    %store normalized trial weights
    for m=1:rogers2024.trialStructure.nSessions
        normTrialWeights(:,a,m) = triFac{a}(:,domFacs(a,m)) ./ max(triFac{a}(:,domFacs(a,m)));
    end
end

%plot
colors = {[.5 .5 .5],'k','r','g','b'};
figure
for m=1:rogers2024.trialStructure.nSessions
    subplot(5, 1, m)
    for g=1:length(rogers2024.groups)
        data=normTrialWeights(:, rogers2024.groups(g).members, m);
        plot(1:rogers2024.trialStructure.nTrialsTotal, mean(data,2), 'color', colors{g}, 'Linewidth', 1)
        hold on
        y_upper = mean(data, 2) + std(data, [], 2) ./ sqrt(size(data, 2));
        y_lower = mean(data, 2) - std(data, [], 2) ./ sqrt(size(data, 2));
        fill( ...
            [1:rogers2024.trialStructure.nTrialsTotal, fliplr(1:rogers2024.trialStructure.nTrialsTotal)], ...
            [y_upper; flipud(y_lower)], ...
            colors{g}, ...
            'FaceAlpha', 0.2, ...
            'EdgeColor', 'none' ...
        );
        hold on
    end
    ylim([0, 1])
    % ylabel('Normalized trial weights')
    % xlabel('Time')
    % title(strcat(rogers2024.sessions(n).name,'-Dominant Component'))
    % if m==5
    %    legend({rogers2024.groups.name})
    % end
end


%% 11. Plot normalized temporal factors. for Supp Fig. 3F-J

for a=1:length(rogers2024.animals)
    %call trial weights from a given animal    
    %store normalized trial weights
    for m=1:rogers2024.trialStructure.nSessions
        normTempWeights(:,a,m) = timFac{a}(:, domFacs(a,m)) ./ max(timFac{a}(:, domFacs(a,m)));
    end
    
end

%plot
colors = {[.5 .5 .5],'k','r','g','b'};
figure
for m=1:rogers2024.trialStructure.nSessions
    subplot(5, 1, m)
    for g=1:length(rogers2024.groups)
        data=normTempWeights(:, rogers2024.groups(g).members, m);
        plot(1:trialLength, mean(data,2), 'color', colors{g}, 'Linewidth', 1)
        hold on
        y_upper = mean(data, 2) + std(data, [], 2) ./ sqrt(size(data, 2));
        y_lower = mean(data, 2) - std(data, [], 2) ./ sqrt(size(data, 2));
        fill([1:trialLength, fliplr(1:trialLength)], [y_upper; flipud(y_lower)], colors{g}, 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        hold on
    end
    ylim([0, 1])
    % ylabel('Normalized trial weights')
    % xlabel('Time')
    % title(strcat(rogers2024.sessions(n).name,'-Dominant Component'))
    % if m==5
    %     legend(rogers2024.groups.name)
    % end
end

%% 12. Assess cumulative distribution of neuron weights. for Fig. 4G,H

% set threshold range
threshold = linspace(0,10,11);

%define storage vectors for fraction of cells included in Acq-, Ext1-, and
%Ext3-dominant ensembles
nComps = 3;
comps = [2 3 5];
cellsInComp = zeros(length(rogers2024.animals), length(threshold),nComps);

for t = 1:length(threshold)
    aCells = [];
    e1Cells = [];
    e3Cells = [];
    
    for n = 1:length(rogers2024.animals)
        
        %call neuron factor weights from animal
        weights = neuFac{n};
        
        for c = 1:nComps
        %calculate fractions of neurons with w > threshold and store
            cellsInComp(n,t,c) = length(find(weights(:,domFacs(n,comps(c)))>threshold(t)))/nCells(n);
        end
        
        %record ensembles at each threshold
        if n == 1
            c = 0;
        else
            c = sum(nCells(1:n-1));
        end
    
        aCells = [aCells; find(weights(:,domFacs(n,2))>threshold(t))+c];
        e1Cells = [e1Cells; find(weights(:,domFacs(n,3))>threshold(t))+c];
        e3Cells = [e3Cells; find(weights(:,domFacs(n,5))>threshold(t))+c];

    end
    
    %collect all the cells in the Acq, Ext1, and Ext3 ensembles at the 
    %predetermined varying thresholds
    collection{t,1} = aCells;
    collection{t,2} = e1Cells;
    collection{t,3} = e3Cells;
end

%plot
colors = {'r','y','c'};
figure
for c=1:nComps
    data = cellsInComp(:,:,c);
    errorbar(threshold, mean(data), std(data)/sqrt(length(rogers2024.animals)-1), colors{c}, 'Linewidth', 3)
    hold on
end
hold on
xline(1)
xlabel('Threshold weight')
ylabel('Fraction of neurons')
legend({'Acq-Dominant','Ext1-Dominant','Ext3-Dominant'})
title('Choosing dominant neurons')


%% 13. Identify Acq-Dominant, Ext1-Dominant, and Ext3-Dominant ensembles and their overlaps for Fig. 5B, 6A

%set threshTest to 0 if reproducing main figs. set threshTest to 1 if
%reproducing Supp Fig 6. note that null thresholds were not determined for
%non-shock mice
threshTest = 0;

%define groups
responders = rogers2024.groups(3).members;
nonresponders = rogers2024.groups(4).members;
rapid = rogers2024.groups(1).members;
slow = rogers2024.groups(2).members;
nonshock = rogers2024.groups(5).members;

if threshTest == 1
    threshes = table2array(readtable('fMeans.csv'));
    threshold = [threshes([15:21 1:14]); [1.55; 1.28; 1.7; 0.975]];
    nAnimals = 21;
end

%create storage vectors for all Acq-, Ext1-, and Ext3-dominant cells
acqDomCells = [];
ext1DomCells = [];
ext3DomCells = [];

%create storage vectors for their overlaps
acqOnly = [];
ext1Only = [];
ext3Only = [];
acqExt1 = [];
acqExt3 = [];
ext1Ext3 = [];
acqExt1Ext3 = [];

%allocate vectors to assign cell indices to for each group
cR = [];
cNR = [];
cSS = [];
cSR = [];
cNS = [];
for a = 1:nAnimals
    if threshTest == 1
        th = threshold(a);
    else
        th=1;
    end
    %number to add to convert within-animal cell index to pooled cell index
    if a == 1
        c = 0;
    else
        c = sum(nCells(1:a-1));
    end
    
    %load cell indices into their group
    if ismember(a,responders)
        cR = [cR 1+c:nCells(a)+c];
    elseif ismember(a,nonresponders)
        cNR = [cNR 1+c:nCells(a)+c];
    elseif ismember(a,rapid)
        cSR = [cSR 1+c:nCells(a)+c];
    elseif ismember(a,slow)
        cSS = [cSS 1+c:nCells(a)+c];
    elseif ismember(a,nonshock)
        cNS = [cNS 1+c:nCells(a)+c];
    end
    
    %call neuron factor weights
    weights = neuFac{a};
    
    %dominant cells are defined by those weights > 1 in the component
    %dominating that session


    aCells = find(weights(:,domFacs(a,2))>th);
    e1Cells = find(weights(:,domFacs(a,3))>th);
    e3Cells = find(weights(:,domFacs(a,5))>th);
    
    weightMat{a,1} = weights(aCells,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,2} = weights(e1Cells,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,3} = weights(e3Cells,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    
    overlap = aCells(ismember(aCells,e1Cells));
    
    %cells dominating all three sessions
    ae1e3 = overlap(ismember(overlap,e3Cells));
    
    %cells dominating all acq, ext1
    ae1 = overlap(~ismember(overlap,e3Cells));
    
    %cells dominating all acq, ext3
    overlap = aCells(ismember(aCells,e3Cells));
    
    ae3 = overlap(~ismember(overlap,e1Cells));
    
    %cells dominating all ext1, ext3
    overlap = e1Cells(ismember(e1Cells,e3Cells));
    mat(a,1) = length(overlap)./nCells(a);
    e1e3 = overlap(~ismember(overlap,aCells));
    
    %cells dominating only acq
    aOnly = aCells(~ismember(aCells,ae1e3));
    aOnly = aOnly(~ismember(aOnly,ae1));
    aOnly = aOnly(~ismember(aOnly,ae3));
    
    %cells dominating only ext1
    e1Only = e1Cells(~ismember(e1Cells,ae1e3));
    e1Only = e1Only(~ismember(e1Only,ae1));
    e1Only = e1Only(~ismember(e1Only,e1e3));
    
    %cells dominating only ext3
    e3Only = e3Cells(~ismember(e3Cells,ae1e3));
    e3Only = e3Only(~ismember(e3Only,e1e3));
     e3Only = e3Only(~ismember(e3Only,ae3));
     
    %store fraction of overlap for individual animals
    ensOverlap(a,1,1) = length(aOnly)/length(aCells);
    ensOverlap(a,2,1) = length(ae1)/length(aCells);
    ensOverlap(a,3,1) = length(ae3)/length(aCells);
    ensOverlap(a,4,1) = length(ae1e3)/length(aCells);
    
    ensOverlap(a,1,2) = length(e1Only)/length(e1Cells);
    ensOverlap(a,2,2) = length(ae1)/length(e1Cells);
    ensOverlap(a,3,2) = length(e1e3)/length(e1Cells);
    ensOverlap(a,4,2) = length(ae1e3)/length(e1Cells);
    
    ensOverlap(a,1,3) = length(e3Only)/length(e3Cells);
    ensOverlap(a,2,3) = length(ae3)/length(e3Cells);
    ensOverlap(a,3,3) = length(e1e3)/length(e3Cells);
    ensOverlap(a,4,3) = length(ae1e3)/length(e3Cells);
    
    %assign neurons indices for activity pooled across animals and store in
    %vectors
    weightMat{a,4} = weights(aOnly,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,5} = weights(e1Only,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,6} = weights(e3Only,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,7} = weights(ae1,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,8} = weights(ae3,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,9} = weights(e1e3,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    weightMat{a,10} = weights(ae1e3,[domFacs(a,1) domFacs(a,2) domFacs(a,3) domFacs(a,4) domFacs(a,5)]);
    
    ensAns(a,1).acqcells = aCells+c;
    ensAns(a,1).ext1cells = e1Cells+c;
    ensAns(a,1).ext3cells = e3Cells+c;
    ensAns(a,1).a = aOnly+c;
    ensAns(a,1).e1 = e1Only+c;
    ensAns(a,1).e3 = e3Only+c;
    ensAns(a,1).ae1 = ae1+c;
    ensAns(a,1).ae3 = ae3+c;
    ensAns(a,1).e1e3 = e1e3+c;
    ensAns(a,1).ae1e3 = ae1e3+c;
    
    acqDomCells = [acqDomCells; aCells + c];
    ext1DomCells = [ext1DomCells; e1Cells + c];
    ext3DomCells = [ext3DomCells; e3Cells + c];
    
    acqOnly = [acqOnly; aOnly+c];
    ext1Only = [ext1Only; e1Only+c];
    ext3Only = [ext3Only; e3Only+c];
    acqExt1 = [acqExt1; ae1+c];
    acqExt3 = [acqExt3; ae3+c];
    ext1Ext3 = [ext1Ext3; e1e3+c];
    acqExt1Ext3 = [acqExt1Ext3; ae1e3+c];
end
% 
% 
labels = {'Overlap with Acq-Dom Neurons','Overlap with Ext1-Dom Neurons','Overlap with Ext3-Dom Neurons'};
labels2{1} = ['Acq Only', 'Acq/Ext1', 'Acq/Ext3', 'Acq/Ext1/Ext3'];
labels2{2} = ['Ext1 Only', 'Acq/Ext1', 'Ext1/Ext3', 'Acq/Ext1/Ext3'];
labels2{3} = ['Ext3 Only', 'Acq/Ext3', 'Ext1/Ext3', 'Acq/Ext1/Ext3'];

nComps=3;
figure
for n = 1:nComps
    for m = 1:length(rogers2024.groups)
        subplot(3,4,m+4*(n-1))
        barWithError(1:4, ensOverlap(rogers2024.groups(m).members,:,n), 1)
        ylabel('Fraction of Cells')
        xticks(1:4)
        xticklabels(labels2{n})
        ylabel(labels{n})
        ylim([0 .8])
        if n==1
        title(rogers2024.groups(m).name)
        end
    end
end

%% 14. Use a Fisher linear decoder to predict group based on average activity (z-score relative to Acquisition) of each ensemble over time for fig. 6B

%TO AVOID HAVING TO RE-RUN PREVIOUS SCRIPTS, LOAD THESE FILES WITH
%ENSEMBLES' CELL INDICES
%
% IF RUNNING THRESHTEST SKIP THIS SECTION
% 
% poolMat = [];
% 
% for n=1:length(rogers2024.animals)
%     poolMat = [poolMat reshape(tensors{n},60*34,nCells(n))];
% end

%define time windows
timeHab = 1:trialLength*max(rogers2024.sessions(1).trials);
timeAcq = trialLength*max(rogers2024.sessions(1).trials)+1:trialLength*max(rogers2024.sessions(2).trials);
timeExt1 = trialLength*max(rogers2024.sessions(2).trials)+1:trialLength*max(rogers2024.sessions(3).trials);
timeExt2 = trialLength*max(rogers2024.sessions(3).trials)+1:trialLength*max(rogers2024.sessions(4).trials);
timeExt3 = trialLength*max(rogers2024.sessions(4).trials)+1:trialLength*max(rogers2024.sessions(5).trials);

%normalize activity to acquisition
poolMat = (poolMat-mean(poolMat(481:960,:)))./std(poolMat(481:960,:));


%define ensembles to loop through
ens = fieldnames(ensAns(1));



%create matrices of average ensembles acts over time in the session of
%interest
sess = 'Ex1'; %set to Acq, Ex1, or Ex3

if strcmp(sess, 'Ex1')
    times = timeExt1;
    s=3;
elseif strcmp(sess, 'Ex3')
    times = timeExt3;
    s=5;
elseif strcmp(sess, 'Ex2')
    times = timeExt2;
    s=4;
elseif strcmp(sess, 'Hab')
    times = timeHab;
    s=1;
elseif strcmp(sess, 'Acq')
    times = timeAcq;
    s=2;
end
timeCell= {timeHab, timeAcq, timeExt1, timeExt2, timeExt3};
%define binary class vector
classes = [zeros(length(times),1); ones(length(times),1)];

ens = fieldnames(ensAns(1));
numEnsembles = length(ens);

sepActs = cell(25,10,2,5);
set1 = zeros(10,360);
set2 = zeros(25,10,360);
set3 = zeros(25,10,360);
for a=1:25
    if a==1
       c=0;
    else
       c=sum(nCells(1:a-1));
    end
    for m=[3 5]
        for l = 1:numEnsembles+1
            if l<11
            pop = ensAns(a).(ens{l});
            
                
            x=0;
            psave = [];
            if ismember(a,rapid)
                if ismember(l,[3,6])
                for p=1:length(pop)
                if mean(poolMat(timeExt3,pop(p)))>8
                    psave = [psave p];
                    disp(psave)
                end
                end
                end
                
                if ~isempty(psave)
                if ~isempty(pop)
                    pop(psave) = [];
                end
                end
            end
            
            
            %pop(ismember(pop,freezeEncoding{a,1}+c)) = [];
            idx=find(fzHz{a,m}==1);
            sepActs{a,l,1,m} = mean(poolMat(timeCell{m}(idx),pop),2);
            idx=find(fzHz{a,m}==0);
            sepActs{a,l,2,m} = mean(poolMat(timeCell{m}(idx),pop),2);
            sepActs{a,l,3,m} = mean(poolMat(timeCell{m},pop),2);
            lens(a,m,l,1) =  size(sepActs{a,l,1,m},1);
            lens(a,m,l,2) =  size(sepActs{a,l,2,m},1);
            lens(a,m,l,3) =  size(sepActs{a,l,3,m},1);
            else
                for ll = 1:10
                    set1(ll,:) = sepActs{a,ll,3,m};
                end
                sepActs{a,l,3,m} = set1;
                lens(a,m,l,1) =  size(sepActs{a,l,1,m},1);
                lens(a,m,l,2) =  size(sepActs{a,l,2,m},1);
                lens(a,m,l,3) =  size(sepActs{a,l,3,m},1);
            end
        end
    end
end

%number of iterations of decoder
numIt = 100;

%set of comparisons - res vs. nonres, res vs. sal, nonres vs. sal
comparisons = 0:2;

%store 100 iterations of the model (rows), for each ensemble plus a round
%where every ensemble is a predictor (columns), for each comparison
%(columns')
%accuracy = zeros(numIt,numEnsembles+1,2);
%accuracyShuff = zeros(numIt,numEnsembles+1,2,3,2);
sOI = [3 5];
for f=3
    for s=1:2
    for l = 1:10
        g1 = rogers2024.groups(3).members;
        g2 = rogers2024.groups(4).members;
        g3 = rogers2024.groups(1).members;
  
        
        if l==11
            ensActs1 = [];
            for a=1:length(g1)
                set = sepActs{g1(a),l,f,sOI(s)};
                ensActs1 = [ensActs1 set];
            end
            
            ensActs2 = [];
            for a=1:length(g2)
                set = sepActs{g2(a),l,f,sOI(s)};
                ensActs2 = [ensActs2 set];
            end
            
            ensActs3 = [];
            for a=1:length(g3)
                set = sepActs{g3(a),l,f,sOI(s)};
                ensActs3 = [ensActs3 set];
            end  
        
            mL = min([length(ensActs1) length(ensActs2) length(ensActs3)]);
  
        
        data = [ensActs1(:,(length(ensActs1)-mL)/2+1:end-(length(ensActs1)-mL)/2)'; ensActs2(:,(length(ensActs2)-mL)/2+1:end-(length(ensActs2)-mL)/2)'; ensActs3(:,(length(ensActs3)-mL)/2+1:end-(length(ensActs3)-mL)/2)']; %compare responders to nonresponders
        data(isnan(data))=0;
        
        else
        
        ensActs1 = [];
        mL = min(lens(g1,sOI(s),l,f));
        for a=1:length(g1)
            set = sepActs{g1(a),l,f,sOI(s)};
            ensActs1 = [ensActs1; set((length(set)-mL)/2+1:end-(length(set)-mL)/2)];
        end
        %ensActs1 = mean(ensActs1)
        ensActs2 = [];
        mL = min(lens(g2,sOI(s),l,f));
        for a=1:length(g2)
            set = sepActs{g2(a),l,f,sOI(s)};
            ensActs2 = [ensActs2; set((length(set)-mL)/2+1:end-(length(set)-mL)/2)];
        end
        
        ensActs3 = [];
        mL = min(lens(g3,sOI(s),l,f));
        for a=1:length(g3)
            set = sepActs{g3(a),l,f,sOI(s)};
            ensActs3 = [ensActs3; set((length(set)-mL)/2+1:end-(length(set)-mL)/2)];
        end
        
        mL = min([length(ensActs1) length(ensActs2) length(ensActs3)]);
  
        
        data = [ensActs1((length(ensActs1)-mL)/2+1:end-(length(ensActs1)-mL)/2); ensActs2((length(ensActs2)-mL)/2+1:end-(length(ensActs2)-mL)/2); ensActs3((length(ensActs3)-mL)/2+1:end-(length(ensActs3)-mL)/2)]; %compare responders to nonresponders
        data(isnan(data))=0;
        
        end
        
        
        lD = length(data)/3;
        for c=1:3
            if c==1
                dat = data(1:2*lD);
                
            elseif c==2
                dat = [data(1:lD); data(1+2*lD:end)];
            else
                dat = data(1+lD:end);
            end
            
            for n=1:100
                classes = [zeros(mL,1); ones(mL,1)];
                [X,Y,x,y] = split_data(dat,classes,.5); %split_data is a custom function, described below that randomly selects test_size % of your data to test on and 1 - test_size to train on
                model = fitcdiscr(x,y); %fit discr fits a linear model to your training data & classes

                predictedLabels  = predict(model, X); %predict applies your model to your test data to generate class predictions

                %count true labels in your test set to normalize confusion
                %matrix
                normVec = zeros(length(unique(y)),1);
                for k=1:length(unique(y))
                    normVec(k,1) = sum(Y==k-1);
                end

                %different sessions have different numbers of unique behaviors 
                accuracy(n,l,f,c,1,s) = mean(predictedLabels == Y);

                shuffle = randperm(length(Y));
                Y = Y(shuffle);
                
                predictedLabels  = predict(model, X); %predict applies your model to your test data to generate class predictions
                accuracy(n,l,f,c,2,s) = mean(predictedLabels == Y); %the accuracy of your model is the number of instances in which your predicted classes matched your a priori classes divided by total number of predictions

            end
        end
    end
    end
end

titles = {'Acq-Dom','Ext1-Dom','Ext3-Dom','Acq Only', 'Ext1 Only', 'Ext3 Only', 'Acq/Ext1', 'Acq/Ext3', 'Ext1/Ext3', 'Acq/Ext1/Ext3','All ensembles'};

tits = {'Decoding between responders and nonresponders during Ext3','Decoding between responders and rapid mice during Ext3','Decoding between nonresponders and rapid mice during Ext3'};
figure; 
for c=1:3
    subplot(3,1,c)
barWithError(1:10, squeeze(accuracy(:,1:10,2,c,1,2)), 0.5)
hold on
barWithError(1:10, squeeze(accuracy(:,1:10,1,c,1,2)), 0.5)
hold on
barWithError(1:10, squeeze(accuracy(:,1:10,3,c,1,2)), 0.5)
hold on
ylim([.5,1])
legend({'Motion','','Freezing'})
xticklabels(titles(1:10))
ylabel('Accuracy')
title(tits{c})
xtickangle(45)
end

tits = {'Decoding between responders and nonresponders during Ext1','Decoding between responders and rapid mice during Ext1','Decoding between nonresponders and rapid mice during Ext1'};
figure; 
for c=1:3
    subplot(3,1,c)
barWithError(1:10, squeeze(accuracy(:,1:10,2,c,1,1)), 0.5)
hold on
barWithError(1:10, squeeze(accuracy(:,1:10,1,c,1,1)), 0.5)
hold on
barWithError(1:10, squeeze(accuracy(:,1:10,3,c,1,1)), 0.5)
hold on
ylim([.5,1])
legend({'Motion','','Freezing'})
xticklabels(titles(1:10))
ylabel('Accuracy')
title(tits{c})
xtickangle(45)
end

%% 15. Extract average activity (z-score relative to Acquisition) of each cell in the ensemble during Extinction 1, Extinction 3. for fig. 5D-F, fig. 6C-I, fig. 7B-D

for d=1:3
for n = 1:25
    if n==1
       c=0;
    else
       c=sum(nCells(1:n-1));
    end
    %loop through all ensembles
    for l=1:10

        cellPop = ensAns(n).(ens{l});
        %cellPop(ismember(cellPop,unique([freezeEncoding{n,1}])+c)) = [];
        %calculate average activity of cells during Extinction 1 and
        %Extinction 3
        %cellPop(ismember(cellPop,freezeEncoding{n,3}+c))= [];
        %(fzHz{n,5}==0)
        %(fzHz{n,3}==0)
         if d==1
             idx1 = (fzHz{n,3}==0);
             idx2 = (fzHz{n,5}==0);
         elseif d==2
             idx1 = (fzHz{n,3}==1);
             idx2 = (fzHz{n,5}==1);
         else
             idx1 = 1:360;
             idx2 = 1:360;
         end
        x =  [mean(poolMat(timeExt1(idx1),cellPop))' mean(poolMat(timeExt3(idx2),cellPop))']; % fzHz{n,5}==0)
        
        
        %store activity matrices
        changedActivity2{l,n,d} = x;
       
    end
    
end
%pool cells

for n=1:4
    group = rogers2024.groups(n).members;
    for l=1:10
        pop=[];
        pgroup=[];
        for g=1:length(group)
            an = group(g);
       
            if l==1
                pop=[pop; changedActivity2{4,an,d}; changedActivity2{7,an,d}; changedActivity2{8,an,d}; changedActivity2{10,an,d}];
            elseif l==2
                pop=[pop; changedActivity2{5,an,d}; changedActivity2{7,an,d}; changedActivity2{9,an,d}; changedActivity2{10,an,d}];
            elseif l==3
                  pop=[pop; changedActivity2{6,an,d}; changedActivity2{5,an,d}; changedActivity2{9,an,d}; changedActivity2{10,an,d}]; %
            else
                pop = [pop; changedActivity2{l,an,d}];
            end
            %pop = [pop; changedActivity2{l,an}];
            pgroup = [pgroup; ensAns(an).(ens{l})];
        end
        
    
    popGroup{l,n} = pgroup;
    pop(isnan(pop)) = 0;
    changedActivity{l,n,d} = pop;    
 
    end

end
sOI=1:10;
titles = {'Acq','Ext1','Ext3','Acq-Only','Ext1-Only','Ext3-Only','Acq/Ext1','Acq/Ext3','Ext1/Ext3','Acq/Ext1/Ext3'};
col = {[.5 .5 .5],'k', 'c','r'};%,
grOI = [1 2 3 4];
figure
for l=1:10
subplot(1,10,l)
for g=1:4
    if size(changedActivity{sOI(l),grOI(g),d},1)<2
        continue
    else
        errorbar(1:2, mean(changedActivity{sOI(l),grOI(g)}), std(changedActivity{sOI(l),grOI(g)})./(sqrt(length(changedActivity{sOI(l),grOI(g)}))-1), 'Linewidth', 3, 'Color', col{g})
    end
hold on
end
xlim([0 3])
xticks(1:2)
xticklabels({'Ext1','Ext3'})
title(titles{sOI(l)})
ylabel('Change in Activity (zscore)')
end
end

for n=1:25
    pop=ensAns(n).acqcells;
    %pop2 = sigCells{n,3,4}{1};
    mA(n) = mean([changedActivity2{4,n}(:,1); changedActivity2{7,n}(:,1); changedActivity2{8,n}(:,1); changedActivity2{10,n}(:,1)]);
    %pop=ensAns(n).ext3cells;
    pop=mean([changedActivity2{6,n}(:,2); changedActivity2{8,n}(:,2); changedActivity2{9,n}(:,2); changedActivity2{10,n}(:,2)]);
    if ismember(n,rogers2024.groups(1).members)
        pop(pop>8) = [];
    end
        
    mE(n) = mean(pop);
end
mA(isnan(mA)) = 0;
mE(isnan(mE)) = 0;



%% 16. TCA reconstructions using representative tca neurons from animal 18
an=18;
weights = neuFac{an};
pop = ensAns(an).a-sum(nCells(1:an-1));
fac = domFacs(an,2);
maxW = max(weights(pop,fac));

neu(1) = find(weights(:,fac)==maxW);

pop = ensAns(an).e1-sum(nCells(1:an-1));
fac = domFacs(an,3);
maxW = max(weights(pop,fac));

neu(2) = find(weights(:,fac)==maxW);

pop = ensAns(an).e3-sum(nCells(1:an-1));
fac = domFacs(an,5);
maxW = max(weights(pop,fac));

neu(3) = find(weights(:,fac)==maxW);

pop = ensAns(an).ae1-sum(nCells(1:an-1));
fac = domFacs(an,2);
maxW = max(weights(pop,fac));

neu(4) = find(weights(:,fac)==maxW);

pop = ensAns(an).ae3-sum(nCells(1:an-1));
fac = domFacs(an,2);
maxW = max(weights(pop,fac));

neu(5) = find(weights(:,fac)==maxW);

pop = ensAns(an).e1e3-sum(nCells(1:an-1));
fac = domFacs(an,3);
maxW = max(weights(pop,fac));

neu(6) = find(weights(:,fac)==maxW);

pop = ensAns(an).ae1e3-sum(nCells(1:an-1));
fac = domFacs(an,2);
maxW = max(weights(pop,fac));

neu(7) = find(weights(:,fac)==maxW);


dat = poolMat(timeAcq,sum(nCells(1:an-1))+1:sum(nCells(1:an)));
datA = zeros(trialLength,nCells(an),8);
for t=1:8
    datA(:,:,t) = dat(1+trialLength*(t-1):trialLength*t,:);
end
datA = zscore(datA);


dat = poolMat(timeExt1,sum(nCells(1:an-1))+1:sum(nCells(1:an)));
datE1 = zeros(trialLength,nCells(an),6);
for t=1:6
    datE1(:,:,t) = dat(1+trialLength*(t-1):trialLength*t,:);
end
datE1 = zscore(datE1);

dat = poolMat(timeExt3,sum(nCells(1:an-1))+1:sum(nCells(1:an)));
datE3 = zeros(trialLength,nCells(an),6);
for t=1:6
    datE3(:,:,t) = dat(1+trialLength*(t-1):trialLength*t,:);
end
datE3 = zscore(datE3);



for n=1:7
    wMat(n,1) = weights(neu(n),domFacs(an,2));
    wMat(n,2) = weights(neu(n),domFacs(an,3));
    wMat(n,3) = weights(neu(n),domFacs(an,5));
end


titles={'AcqOnly','Ext1Only','Ext3Only','Acq/Ext1','Acq/Ext3','Ext1/Ext3','Acq/Ext1/Ext3'};
sOI = [2 3 5];
for m=1:3
tim = rogers2024.sessions(sOI(m)).trials;
figure
for n=1:7
fH = timFac{an}(:,domFacs(an,1)).*weights(neu(n),domFacs(an,1)).*mean(triFac{an}(tim,domFacs(an,1)));
fA = timFac{an}(:,domFacs(an,2)).*weights(neu(n),domFacs(an,2)).*mean(triFac{an}(tim,domFacs(an,2)));
fE1 = timFac{an}(:,domFacs(an,3)).*weights(neu(n),domFacs(an,3)).*mean(triFac{an}(tim,domFacs(an,3)));
fE2 = timFac{an}(:,domFacs(an,4)).*weights(neu(n),domFacs(an,4)).*mean(triFac{an}(tim,domFacs(an,4)));
fE3 = timFac{an}(:,domFacs(an,5)).*weights(neu(n),domFacs(an,5)).*mean(triFac{an}(tim,domFacs(an,5)));

subplot(2,4,n)
dat=sum([fH fA fE1 fE2 fE3],2);
dat = dat-min(dat);
dat = dat./max(dat);
plot(dat)
hold on
dat1 = traces{n,m}(1:60)-min(traces{n,3}(1:60));
dat1 = dat1./max(dat1);
plot(dat1)
ylim([-.1 1.1])
legend('Reconstructed','Real')
title(titles{n})

[rRec(n,m), pRec(n,m)] = corr(dat,dat1);
pRec(pRec>.05) = 0;
%rRec(pRec==0) = 0;
end
end
rRec=abs(rRec);
figure
subplot(121)
heatmap(rRec,'colorlimits',[0 1])
colormap('parula')
subplot(122)
heatmap(pRec,'colorlimits',[0 .05])
colormap('parula')

%% 17. identify overlap of TCA ensembles with stim-responsive neurons (fig 8, supp 10)
for a=1:25
    if a==1
        c=0;
    else
        c=sum(nCells(1:a-1));
    end
    for l=1:3
        pop1 = ensAns(a).(ens{l})-c;
        for k=1:3
            for m=1:5
                if k<3
                pop2 = [sigCells{a,m,k}{1} sigCells{a,m,k}{2}];
                else
                    pop2 = [sigCells{a,m,4}{1} sigCells{a,m,4}{2}];
                end
                fracEns(a,l,k,m) = mean(ismember(pop1,pop2));
            end
        end
    end
end