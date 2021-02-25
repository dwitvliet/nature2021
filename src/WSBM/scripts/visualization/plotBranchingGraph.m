function plotBranchingGraph(labels, showNodeNames, outputDir)
% given a set of labels, plot out their intersections as branching plot
% labels: n x numDS matrix, showNodenames = 1 if showing node names on each
% group
if nargin == 1
    showNodeNames = 0;
    outputDir = 0;
elseif nargin == 2
    outputDir = 0;
end
%% extract some meta information from the labels
[~, numDS] = size(labels);
info = nan(numDS, 3); % # of K's, start index, end index
i = 0;
for idx = 1:numDS
    numK = numel(unique(labels(:, idx)));
    info(idx,:) = [numK i+1 i+numK];
    i = i+numK;
end
numNodes = info(size(info, 1), size(info, 2));
%% create community overlap
overlaps = cell(1, numDS -1);
for t = 1:numel(overlaps)
    m = info(t, 1);
    n = info(t+1, 1);
    overlaps{t} = nan(m, n);
    for i = 1:m
        for j = 1:n
            same = labels(:, t) == i & labels(:, t+1) == j;
            overlaps{t}(i,j) = sum(same(:));
        end
    end
end

%% adjacency matrix
adjMat = zeros(info(size(info, 1), size(info, 2)));
for t = 1:numel(overlaps)
    adjMat(info(t, 2):info(t, 3), info(t+1, 2):info(t+1, 3)) = overlaps{t};
end


[~, parents] = max(adjMat,[], 1);
parents(info(1, 2):info(1, 3)) = 0;


% rearrange adjacency matrix so that the children are put in a position
% corresponds to the parent.
numIter = 0;
while ~all(parents == sort(parents))
    [~, idx] = sort(parents);
    adjMat = adjMat(idx, idx);
    
    [~, parents] = max(adjMat,[], 1);
    parents(info(1, 2):info(1, 3)) = 0;
    
    numIter = numIter + 1;
    if numIter > 300
        disp("too many iterations and it hasn't converge");
        break
    end
end

% manual override
idx = [2 1 4 3 5 7 6 8 10 9 11 12 14 13 15 16 19 17 18 20 21 22 25 23 24 26 27 28 31 29 30 32 33 34];
adjMat = adjMat(idx, idx);
colors = [ ...
    [0.9333 0.9333 0.9333]; [0.9333 0.9333 0.9333]; ...
    [1.0000 0.6824 0.6667]; [0.8118 0.4941 1.0000]; [0.5137 0.6941 1.0000]; ...
    [1.0000 0.6824 0.6667]; [0.8118 0.4941 1.0000]; [0.5137 0.6941 1.0000]; ...
    [1.0000 0.6824 0.6667]; [0.8118 0.4941 1.0000]; [0.5137 0.6941 1.0000]; [0.3843 0.6118 1.0000]; ...
    [1.0000 0.6824 0.6667]; [0.8118 0.4941 1.0000]; [0.5137 0.6941 1.0000]; [0.3843 0.6118 1.0000]; ...
    [1.0000 0.6824 0.6667]; [0.8824 0.6863 1.0000]; [0.8118 0.4941 1.0000]; [0.6980 0.8118 1.0000]; [0.5137 0.6941 1.0000]; [0.3843 0.6118 1.0000]; ...
    [1.0000 0.6824 0.6667]; [0.8824 0.6863 1.0000]; [0.8118 0.4941 1.0000]; [0.6980 0.8118 1.0000]; [0.5137 0.6941 1.0000]; [0.3843 0.6118 1.0000]; ...
    [1.0000 0.6824 0.6667]; [0.8824 0.6863 1.0000]; [0.8118 0.4941 1.0000]; [0.6980 0.8118 1.0000]; [0.5137 0.6941 1.0000]; [0.3843 0.6118 1.0000]; ...
];



%% percentage that a community preseve it's parent's composition
out = sum(adjMat, 2);
out = repmat(out, 1, numNodes);
percentage = adjMat * 100 ./ out;
percentage(isnan(percentage)) = 0;



%% visualize the change of communities


%%% calculate x and y position of each node
figure();
x = zeros(1, numNodes);
y = zeros(1, numNodes);
% nodeLabels_disp = cell(size(nodeLabels));
for i = 1:numDS
    x(info(i, 2): info(i, 3)) = i;
    pts = linspace(0, max(info(:, 1)), info(i, 1) + 2);
    y(info(i, 2): info(i, 3)) = pts(2:end-1);
end

G = digraph(adjMat);
h = plot(G, 'XData',x,'YData',y);
markerSize = sum(adjMat, 2);
markerSize(markerSize == 0) = sum(adjMat(:, info(numDS, 2):info(numDS, 3)),1);
%disp(markerSize)
h.MarkerSize = markerSize * 0.5;
h.NodeLabel = [];
h.NodeColor = colors;
%h.NodeLabel = markerSize;
%h.EdgeLabel = nonzeros(adjMat');
h.LineWidth = nonzeros(adjMat')*0.5;
h.EdgeCData = nonzeros(percentage');
h.EdgeAlpha = 1;
h.ArrowSize = 0;
%h.ArrowSize = rescale(nonzeros(adjMat'), 10, 20);



xlabel("Datasets");
ax = gca;
% ylim(ax, [0.3 5.3]);
colormap(flipud(gray))
c = colorbar;
set(c,'position',[.85 .3 .02 .4])

pos = get(c,'Position');
c.Label.String = compose("Percentage\ninherited");
c.Label.Position = [pos(1)/2 pos(2)+115]; 
c.Label.Rotation = 0; 
caxis([0 100])

set(gca, 'XColor', 'black', 'YColor', 'white', 'TickDir', 'out');
set(gca,'box','off')

if outputDir
    fprintf('saving results in %s\n', outputDir);
    savefig(sprintf("%s/analysis_result/plots/cluster_change/%s.fig", result_dir, filename)); 
end
end

