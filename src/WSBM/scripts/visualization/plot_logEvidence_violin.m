function [logEvidences, maxLogEvi, maxLogEviIdx] = plot_logEvidence_violin(statDir, graphDir, datasetLabels)
%PLOT_LOGEVIDENCE_VIOLIN Summary of this function goes here
%   Detailed explanation goes here
logEvidences = load(fullfile(statDir, "logEvidences.mat"));
logEvidences = logEvidences.all_logEvidences;

numDS = size(logEvidences, 1);
fig1=figure();
fig1.Renderer='Painters';
hold on
maxLogEvi = cell(numDS, 3);
maxLogEviIdx = cell(numDS, 3);
DS_k = cell(numDS, 2);
for i = 1:numDS
    if isstruct(logEvidences{i, 1})
        fileNames = {logEvidences{i, 1}.name};
    else
        fileNames = logEvidences{i, 1};
    end
    [fileNames, fileNameIdx] = natsort(fileNames);
    subplot(ceil(numDS/2), 2, i);
    dsIdx = unique(regexp([fileNames{:}], "(?<=_ds)\d(?=_k)", "match"));
    DS_k{i, 1} = str2num(dsIdx{:});
    ks =regexp([fileNames{:}], "(?<=_k)\d\d?(?=.mat)", "match");
    DS_k{i, 2} = ks;
    logEvidence = logEvidences{i, 2};
    logEvidence = logEvidence(fileNameIdx, :);
    violinplot(logEvidence', ks, 'ViolinAlpha', 1, 'ShowData', false);
    [y, x] = max(mean(logEvidence, 2));
    maxLogEvi{i, 3} = y;
    maxLogEviIdx{i, 3} = str2num(ks{x});
    title(sprintf("%s", datasetLabels{i}));
    xlabel("Number of modules");
    ylabel("Log-likelihood score");
end
hold off;

maxLogEvi(:, 1:2) = DS_k;
maxLogEviIdx(:, 1:2) = DS_k;
set(gcf, 'Position',  [100, 100, 700, 1200]);
savefig(sprintf("%s/logEvidence.fig", graphDir));
saveas(gcf, sprintf("%s/logEvidence.pdf", graphDir));
saveas(gcf, sprintf("%s/logEvidence.png", graphDir));
end

