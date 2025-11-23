

function [longRegistered, coordinates] = crossSeshAuto(cellreg, calciums, cellprops)
%this function aligns longitudinally registered datasets
%INPUTS:    cellreg - the cell registration file from IDPS
%           calciums - a cell containing the cell trace files from each session in rows
%OUTPUTS:   longRegistered = a cell of matrices of longitudinally registered cell
%               traces from each session
%
%
%
%Written by Sophie A. Rogers, Corder Laboratory, University of Pennsylvania

    vec = unique(cellreg.global_cell_index);
    cellIDs = [];
    for n = 1:length(vec)
        ind = find(cellreg.global_cell_index==vec(n));
        if length(ind)==length(calciums)
            cellIDs = [cellIDs; cellreg.local_cell_index(ind)'];
        end
    end
   
    for n=1:length(calciums)
        if size(calciums{n}, 2) - 1 < cellIDs(end, n)
            cellIDs = cellIDs(1:end-1, :);
        end
    end
    
    
    
    %extract the traces of registered cells from each session
    for n=1:length(calciums)
        a = unique(cellIDs(:,n))+1;
        f = table2array(calciums{n});
        f = f(:,2:end);
        f(isnan(f)) = 0;
        longRegistered{n} = f(:,a);
    end

    if ~isempty(cellprops)
        coordinates = [cellprops.CentroidX(a) cellprops.CentroidY(a) cellprops.Size(a)];
    elseif isempty(cellprops)
        coordinates = [];
    end
    
    disp(size(longRegistered{1}))
end
