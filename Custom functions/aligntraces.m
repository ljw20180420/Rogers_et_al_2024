function allAligned = aligntraces(traces, stimulusTimes, nSessions, samplingRate, trialStart, trialEnd)
%this function aligns recorded neural activity to stimuli
%INPUTS =   traces - a cell with neural traces from each session
%           stimulusTimes - a cell with stimulus times from each session
%           nSessions - number of sessions
%           samplingRate - sampling rate of video
%           trialStart - amount of time prior to stimulus start to include
%           trialEnd - amount of time after stimulus start to include
%OUTPUTS =  alignedTraces - a cell with matrices of aligned neural traces
%           concatenated in columns
%
%
%
%Written by Sophie A. Rogers, Corder Laboratory, University of Pennsylvania
%%  
    for m=1:nSessions        
        acts = zeros(round(size(traces{m}, 1) / samplingRate), size(traces{m}, 2));
        for ii = 1:round(size(traces{m}, 1) / samplingRate) - 1
            acts(ii,:) = sum(traces{m}((ii-1) * samplingRate + 1:ii * samplingRate, :)) / samplingRate;
        end
        
        for t=1:length(stimulusTimes{m})
            allAligned{t, m} = acts(round(stimulusTimes{m}(t)) - trialStart:round(stimulusTimes{m}(t)) + trialEnd, :);
        end        
    end
    allAligned = vertcat(allAligned{:});
end
