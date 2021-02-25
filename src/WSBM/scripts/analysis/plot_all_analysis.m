% given a folder with model and labels (unique for each ds and k)
% plot some useful analysis results, including
% 1. logEvidence of each data point
% 2. # of pairs for each data point
% 3. branching diagram showing how consecutive results relates to each
% other
% 4. violin plot for mean KS against empirical distribution
run('setup.m');
% inputDir = './results/A-all-postemb';
inputDir = './results/B-all';
% inputDir = './results/C-no-noise';
% inputDir = './results/D-only_stable';
% inputDir = './results/E-all-pruned-to-l1';
% inputDir = './results/F-all-grown-to-adult-random';
% inputDir = './results/G-all-grown-to-adult';
consensusDir = sprintf('%s/consensus', inputDir); % dir of the model to visualize
graphDir = sprintf("%s/analysis_results", inputDir);
statDir = sprintf('%s/stats', inputDir);
dataDir = './data'; %dir storing the empirical dist
datasetLabels = {'L1 (0 hrs)','L1 (5 hrs)','L1 (8 hrs)','L1 (16 hrs)','L2 (23 hrs)','L3 (27 hrs)', 'Adult (50 hrs)', 'Adult (50 hrs)'};
mkdir(graphDir);
% plot log evidence for each datasets
[all_logEvidence, maxLogEvi, maxLogEviIdx] = plot_logEvidence_violin(statDir, graphDir, datasetLabels);
k = [maxLogEviIdx{:, 3}] % the chosen k for each dataset
%% brancing plot
DS = [maxLogEviIdx{:, 1}]; % the datasets to show -> can change to the combination you wanna visualize

k=[2 3 3 4 4 6 6 6];
disp('the datasets to show are');
disp(DS);
disp("the k's to show are");
disp(k);
%%% load in labels and models
numDS = length(DS);
if numDS ~= length(k)
    error('DS and k not equal in length');
end
models = cell(numDS, 1);
labels = nan(212, numDS);

for i = 1:numDS
    model_fileName = sprintf('%s/consensus_model_ds%d_k%d.mat', consensusDir, DS(i), k(i));
    label_fileName = sprintf('%s/consensus_label_ds%d_k%d.mat', consensusDir, DS(i), k(i));

    model = load(model_fileName);
    label = load(label_fileName);
    model = getTheOnlyField(model);
    label = getTheOnlyField(label);
    models{i} = model;
    labels(:, i) = label;

end
%% write the optimal grouping to file
fileID = fopen(sprintf("%s/groupings.csv", graphDir),'w');
for i = 1:numDS
    fprintf(fileID, 'ds = %d,k= %d\n',DS(i), k(i));
    groupLabels = getGroupNodeLabels(viewres(labels(:, i)), 1);
    disp(groupLabels)
    for j = 1: k(i)
        fprintf(fileID, "%s\n", groupLabels{j});
    end
end
fclose(fileID);
%%
plotBranchingGraph(labels, 1);
set(gcf, 'Position',  [100, 100, 1000, 700])
set(gca,'ytick',[])
set(gca,'yticklabel',[])
xticks(1:8)
xticklabels(datasetLabels);
xtickangle(45);
savefig(sprintf("%s/branchingPlot.fig", graphDir));
saveas(gcf, sprintf("%s/branchingPlot.pdf", graphDir));
saveas(gcf, sprintf("%s/branchingPlot.png", graphDir));

