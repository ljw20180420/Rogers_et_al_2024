function changedCells = isSig(traces,baselineTime,targetTime)
%this function checks whether cells are significantly modulated by certain
%   stimuli
% INPUT:    traces  - cell traces over trial of interest
%           baselineTime - array of times considered baseline
%           targetTime - array of times of interest
% OUTPUT:   changedCells - cell containing cells that are up, down, and not
%           regulated during target time
%
%
%
%Written by Sophie A. Rogers, Corder Laboratory, University of Pennsylvania
%%
up = [];
down = [];
non = [];
for n=1:size(traces, 2)
    [P,H] = ranksum(traces(baselineTime,n),traces(targetTime,n));
    if H == 1
        if median(traces(targetTime,n)) > median(traces(baselineTime,n))
            up = [up n];
        else
            down = [down n];
        end
    else
        non = [non n];
    end
end

changedCells{1,1} = up;
changedCells{1,2} = down;
changedCells{1,3} = non;

end
