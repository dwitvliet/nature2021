function [cnsLabel, cnsModel] = fitConsensusModel(resultDir, graphIdx, k, numThread)
    model_path = sprintf("%s/ds%d/models_ds%d_k%d.mat", resultDir, graphIdx, graphIdx, k);
    models = load(model_path);
    models = models.models;
    repeats = numel(models);
    LOOPER_ITER = repeats;
    
    pause(1+rand()); % wtf
    c = parcluster();
    r = rand;
    mkdir(sprintf("/tmp/local_cluster_jobs/ds%d_k%d_%d",graphIdx, k, r));
    c.JobStorageLocation = sprintf("/tmp/local_cluster_jobs/ds%d_k%d_%d",graphIdx, k, r);
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
    
    kLoopBestModels = cell(1, repeats);
    for i = 1: repeats
        kLoopBestModels{i} = models(i);
    end
    kLooperModels = kLoopBestModels;
    
    [ kCentralModel , kSimVec , kSimVUpTri ] = central_model(kLoopBestModels) ; 

    %%% sort the kLooperModels...
    % by aligning all to the 'most central' with Hungarian algo

    % the ref to align to will be the kCentralModel
    [~,ref] = community_assign(kCentralModel) ;

    % preallocate the mat of the bestModels community assignments
    kBestModels_ca_algn = zeros([ kCentralModel.Data.n LOOPER_ITER ]) ;

    % first stack all plausible parcellations
    for idx=1:LOOPER_ITER

        [~,tmp] = ...
            community_assign(kLoopBestModels{idx}.Para.mu) ;

        kBestModels_ca_algn(:,idx) = CBIG_HungarianClusterMatch(ref,tmp) ;

    end

    %%% get kSimVec and LogEvid vec sorted

    % will be size: numKTested x LOOPER_ITER
    kLooperLogEvids = cellfun(@(a) a.Para.LogEvidence,kLooperModels(:,:)) ;
    % get the row with the best evids, and transpose to make col
    kBestLogEvids = kLooperLogEvids' ;

    % sort evid vector and kSim vector
    [ sort_logEvid, sort_logEvidIdx ] = sort(kBestLogEvids,'descend');
    [ sort_kSimVec, sort_kSimVecIdx ] = sort(kSimVec,'descend');

    % show that there is corrleation between logEvid and VI disatnce
    logEvid_kSim_corr = corr(kBestLogEvids,kSimVec)

    %%% lets iterate over the 100 k best models to create consensus model
    % first we have to turn the hungarian reordered ca_align mat into a wsbm
    % prior to be input into wsbm inference fit

    % lets cuttoff the worst 5% of the runs to clear up the prior a little...
    % TODO... maybe no need for this...
    % if no cutting of... then sorting doesn't matter-but Hungarian match still
    % does matter
    onlyKeep = 1%0.95 ;

    % first figure out the cuttoff
    cutoff = floor(LOOPER_ITER * onlyKeep) ; 
    kiter_prior = zeros([ kCentralModel.R_Struct.k kCentralModel.Data.n ]) ;

    % sort the ca_align by similarity
    kBestModels_ca_algn_sort = kBestModels_ca_algn(:,sort_kSimVecIdx);

    % iterate over each node
    for idx=1:(kCentralModel.Data.n)

        kiter_prior(:,idx) = ...
            sum(bsxfun(@eq,kBestModels_ca_algn_sort(idx,1:cutoff), ...
            [1:(kCentralModel.R_Struct.k)]'),2) ./ cutoff ;

    end    
    disp(kiter_prior);
    %%% loop for consensus

    CONSENSUS_ITER = 10 ;
    C = zeros([CONSENSUS_ITER 1]) ;
    kBest = k;
    for idx=1:CONSENSUS_ITER

        runNTimes = 100 ;

        rr = sym_RStruct(kBest) ;
        modIn = { ... 
            'W_Distr', 'Exponential', ...
            'E_Distr', 'Bernoulli', ...
            'alpha', 0.5, ...
            'mainMaxIter', 300 , ...
            'muMaxIter' , 300,...
            'mainTol',0.001, ...
            'muTol', 0.001 ,...
            'verbosity', 0, ...
            'numTrials', 50 ,...
            'mu_0', kiter_prior(:,:,idx)};

        % function [ allModels ] = wsbmFitNTimes( adjMat, rStruct , modelInputs , numFits , numCores)
        cnsnsusModels = wsbmFitNTimes(kCentralModel.Data.Raw_Data,...
            rr,...
            modIn,...
            runNTimes, 16) ;

        tmpCnsnsusCa = zeros([ kCentralModel.Data.n runNTimes ]) ;

        for jdx=1:runNTimes
            [~,tmpCnsnsusCa(:,jdx)] = community_assign(cnsnsusModels(jdx).Model) ;
        end

        tmpAgreeMat = calcConcensusMat(tmpCnsnsusCa);

        % get the consensus consitency
        C(idx) = consensus_consistency(tmpAgreeMat) 

        % make new kiter_prior for new loop
        for kdx=1:(kCentralModel.Data.n)

            kiter_prior(:,kdx,idx+1) = ...
                sum(bsxfun(@eq,tmpCnsnsusCa(kdx,:), ...
                [1:(kCentralModel.R_Struct.k)]'),2) ./ cutoff ;

        end    

        % have we converged? or are we at the end of loop?
        if C(idx) >= 0.95 || idx == CONSENSUS_ITER
            if idx == CONSENSUS_ITER
                warning("didn't reach consensus in %d iterations", CONSENSUS_ITER);
            end
            consensus_kiter_prior = kiter_prior(:,:,idx+1) ;    
            [~,consensus_ca] = community_assign(consensus_kiter_prior) ; 
            consensus_kCentralModel = central_model(cnsnsusModels) ;
            % also add the data back into consensus_kCentral...
            consensus_kCentralModel.Data = kCentralModel.Data ;
            break 
        end

    end

    %%% get centroid of best k and define it as templateModel

    cnsModel = consensus_kCentralModel ; 
    cnsLabel = consensus_ca;
    mkdir(sprintf("%s/consensus/", resultDir));
    [~, txtLabs] = viewres(cnsLabel);
    fid = fopen(sprintf("%s/consensus/txt_consensus_label-ds%d_k%d.txt", resultDir, graphIdx, k),'w');
    fprintf(fid, "Condition: %s, Dataset: %d, k: %d \n", resultDir, graphIdx, k);
    for i = 1:numel(txtLabs)
      fprintf(fid, "Group: %d\n", i);
      fprintf(fid, "%s\n", txtLabs{i});
    end
    fclose(fid);

    % save the inital best model
    save(sprintf("%s/consensus/consensus_model_ds%d_k%d.mat", resultDir, graphIdx, k), "cnsModel");
    save(sprintf("%s/consensus/consensus_label_ds%d_k%d.mat", resultDir, graphIdx, k), "cnsLabel");
end
