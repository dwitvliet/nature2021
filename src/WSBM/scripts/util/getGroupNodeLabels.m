function [groupNodeLabels, sizes] = getGroupNodeLabels(res, oneLine)
    if nargin == 1
        oneLine = 0;
    end
    G = findgroups(res(:, 3));
    sizes = splitapply(@size, G, G);
    sizes = sizes(:, 1);
    groupNodeLabels = splitapply(@aggregateNodes, res(:, 1), G);
    if oneLine
        groupNodeLabels = splitapply(@aggregateNodes_oneLine, res(:, 1), G);
    end
end

function s = aggregateNodes(vec)
% create a string of node with each line consists of at most 5 nodes
    maxNodesPerLine = 5;
    len = size(vec, 1);
    groupings = zeros(len, 1);
    for i = 1: maxNodesPerLine: len
        groupings(i:end, 1) = groupings(i:end, 1) + 1;
    end 
    func1 = @(s) join(s, ', ');
    vec_group = splitapply(func1, vec(:, 1), findgroups(groupings));
    s = join(vec_group, newline);
end

function s = aggregateNodes_oneLine(vec)
% doesn't do the split line business, just concatenate all neuron names to
% one line 
    maxNodesPerLine = 5;
    len = size(vec, 1);
    groupings = zeros(len, 1);
    for i = 1: maxNodesPerLine: len
        groupings(i:end, 1) = groupings(i:end, 1) + 1;
    end 
    func1 = @(s) join(s, ', ');
    vec_group = splitapply(func1, vec(:, 1), findgroups(groupings));
    s = join(vec_group, ',');
end