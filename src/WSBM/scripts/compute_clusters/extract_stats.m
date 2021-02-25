function extract_stats(resultDir)
    % extract information from multiple model (per ds one condition) fitting runs and fit a consensus
    % model
    % Precondition:
    % the input directory is assumed to have the following strcture
    % ./
    %   ds0/
    %       models_ds0_k*.mat
    %       labels_ds0_k*.mat
    %   ds1/
    %       ...
    % stats collected includes
    % 1. LogEvidence for each run
    % 2. number of pairs for each run
    % 3. the consensus matrix for each k for each ds
    % the stats collected will be stored under ./stats/
    
    % set up the read and write directories
    dsDirs = dir(sprintf("%s/ds*", resultDir));
    outputDir = sprintf("%s/stats/", resultDir);
    
    [~, idx]= natsort({dsDirs.name});
    dsDirs = dsDirs(idx);
    % set up the data structure to store the statistics
    fprintf("start to process %s\n", resultDir);
    if numel(dsDirs) == 0
      error("no datasets found");
    else
      fprintf("a total of %d datasets found\n", numel(dsDirs));
    end

    all_logEvidences = cell(numel(dsDirs), 2);
    all_pairs = cell(numel(dsDirs), 2);
%     all_cnsMats = cell(numel(dsDirs), 2);
    for i = 1:numel(dsDirs)
        
        
        dsDir = strcat(resultDir, '/', dsDirs(i).name);
        model_paths = dir(sprintf("%s/models*.mat", dsDir));
        label_paths = dir(sprintf("%s/labels*.mat", dsDir));
        model_paths = natsort({model_paths.name});
        label_paths = natsort({label_paths.name});
        all_logEvidences{i, 1} = model_paths;
        all_pairs{i, 1} = model_paths;
%         all_cnsMats{i, 1} = model_paths;
        
        if numel(model_paths) ~= numel(label_paths)
            error("in ds folder %s, model has %d elements while label has %d elements", dsDir, numel(model_paths), numel(label_paths));
        end
        
        numDS = numel(model_paths);
        logEvidences = nan(numDS, 1);
        pairs = nan(numDS, 1);
%         cnsMats = nan(numDS, 1);

        fprintf("number of models to process in %s is %d\n", dsDirs(i).name, numel(model_paths));
        for k = 1:8 %numel(model_paths)
            model_path = strcat(dsDir, '/', model_paths{k});
            label_path = strcat(dsDir, '/', label_paths{k});
            models = load(model_path);
            models = models.models;
            labels = load(label_path);
            labels = labels.labels;
            repeats = numel(models);
            
            logEvidence = nan(1, repeats);
            pair = nan(1, repeats);
            parfor r = 1:repeats
%                 n = models(r).Data.n;
                logEvidence(r) = models(r).Para.LogEvidence;
                res = viewres(labels(:, r));
                p = countPairs(res);
                pair(r) = p;
            end
%             n = models(1).Data.n;
            logEvidences(k, 1:repeats) = logEvidence;
            pairs(k, 1:repeats) = pair;
%             cnsMats(k, 1:n, 1:n) = calcConcensusMat(labels);
        end
        
        all_logEvidences{i, 2} = logEvidences;
        all_pairs{i, 2} = pairs;
%         all_cnsMats{i, 2} = cnsMats;
    end
    mkdir(outputDir);

    disp("saving outputs");

    save(sprintf("%s/logEvidences.mat", outputDir), 'all_logEvidences');
    save(sprintf("%s/pairs.mat", outputDir),  'all_pairs');
%     save(sprintf("%s/cnsMats.mat", outputDir), 'all_cnsMats');
end
