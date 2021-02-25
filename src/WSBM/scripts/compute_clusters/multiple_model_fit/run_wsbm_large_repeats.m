function [labels, models, LogEvidence, pairs] = run_wsbm_large_repeats(E, k, repeats, outputDir, graphIdx, muIter, mainIter, numThread)
% run wsbm 
    if nargin == 6
      numThread = 1;
    end
    if nargin < 6
      error("not enough input arguments.");
    end

    model_name= sprintf("%s/models_ds%d_k%d.mat", outputDir, graphIdx, k);
    labels_name = sprintf("%s/labels_ds%d_k%d.mat", outputDir, graphIdx, k);


    a = nan(size(E, 1), repeats);
    if exist("models", "var")

        clear("models");
    end

    %pause(1+60*rand()); %wtf
    pause(1+rand());       %wtf
    c = parcluster();
    r = rand;
    mkdir(sprintf("/tmp/local_cluster_jobs/ds%d_k%d_%d",graphIdx, k, r));
    c.JobStorageLocation = sprintf("/tmp/local_cluster_jobs/ds%d_k%d_%d", graphIdx, k, r);
    c.NumWorkers = numThread - 1;
    if exist('parpool')
      p = gcp('nocreate');
    else
      p = parpool(c, numThread - 1);
    end

    if isempty(p)
        poolsize = 0;
    else
        poolsize = p.NumWorkers;
    end
    fprintf("matlab running ds %d k%d with %d workers", graphIdx, k, poolsize);
    tic
    parfor(r = 1: repeats, numThread - 1)
        [label, model] = wsbm(E, k, 'W_Distr', 'Exponential', 'muMaxIter', muIter, 'mainMaxIter', mainIter, 'numTrials', 300);
        a(:, r) = label;
        models(r) = model;
        disp(r);
    end
    labels = a;
    save(model_name, 'models');
    save(labels_name, 'labels');
    toc
end
