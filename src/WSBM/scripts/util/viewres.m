function [res, groupNodeLabels, groupSizes] = viewres(labels, showRes)
    if nargin < 2
        showRes = 0;
    
    end
    map = readtable('./data/node_int.csv');
    map = sortrows(map(1:numel(labels),:), 2);

    groups = addvars(map, labels);
    res = sortrows(groups, [3 1]);
    [groupNodeLabels, groupSizes] = getGroupNodeLabels(res);
    if showRes
        disp(groupNodeLabels);
    end
end

