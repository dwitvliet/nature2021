function [concensusMat, res] = calcConcensusMat(labels, showResults)
    % compute the concensus matrix of a set of clusterings
    % labels is a matrix of size node_num (N) x trial_num(t)
    if nargin == 1
        showResults = 0;
    end
    [N t] = size(labels);
    concensus_mats = nan(N, N, t);
    for i = 1: t
        label = labels(:, i);
        concensus_mats(:,:, i) = (repmat(label, 1, N) == repmat(label', N, 1));
    end
    concensusMat = sum(concensus_mats, 3)/t;
    cMat = 1 - concensusMat;
    method = "weighted";
    tree = linkage(cMat, method);
    k = max(labels(:));
    clustered_labels = cluster(tree, "Maxclust", k);
    res = viewres(clustered_labels);
    if showResults
        reorderedLabels = res(:, 1);
        index = table2array(res(:, 2));
        names = table2array(reorderedLabels);
        figure();
        heatmap(concensusMat(index, index), 'GridVisible', 'off', 'XDisplayLabels', names, 'YDisplayLabels', names)
        title({"Consensus Matrix", sprintf("method: %s, number of clusters: %d, number of pairs: %d", method, k, countPairs(res))});
        colormap jet
    end
end