function [ loopResults , allFitModels, allFitLogEvidence ] = wsbm_looper_wrapper( inputData , modelInputs , loopIters , scoreFunc, numTrialPttrn, priorWeightPttrn, outputDir)

% wrapper for the looper script already in the wsbm code

% initialize the variables we will write into
% so it parallels easily
numModels = numel(modelInputs) ;
loopResults = zeros(numModels, loopIters + 1) ; 
loopResults(1:numModels,1) = 1:numModels ;
%loopResults(numModels + 1, 1) = 999 ;

allFitModels = cell(numModels,loopIters) ;
allFitLogEvidence = nan(loopIters, numModels, numel(numTrialPttrn));

if ~exist('scoreFunc','var') || isempty(scoreFunc)
  scoreFunc = @(model) model.Para.LogEvidence ;
  disp('using log evidence as score func')
else
  disp('using provided score func')
end

if ~exist('numTrialPttrn','var') || isempty(numTrialPttrn)
    numTrialPttrn = [] ;
end

if ~exist('priorWeightPttrn','var') || isempty(priorWeightPttrn)
    priorWeightPttrn = [] ;
end

% parallel_pool = gcp ; 
delete(gcp('nocreate'))

%  pause(1+60*rand());
c = parcluster();
c.NumWorkers = 8;
parpool(c);
% r = rand;
% mkdir(sprintf("%s/local_cluster_jobs/ds%d_k%d_%d", getenv("SCRATCH"),graphIdx, k, r));
% c.JobStorageLocation = sprintf("%s/local_cluster_jobs/ds%d_k%d_%d", getenv("SCRATCH"),graphIdx, k, r);
% c.NumWorkers = numThread;
% if exist('parpool')
%   p = parpool(c);
% else
%   p = parpool(c, numThread);
% end
% 
% if isempty(p)
%     poolsize = 0;
% else
%     poolsize = p.NumWorkers;
% end
ppm1 = ParforProgMon('looper',loopIters,1) ;
parfor idx = 1:loopIters

    %disp('iteration:')
    %disp(idx)
    
    % Fit
    % function [Best_Models,Scores,Models] = wsbmLooper_2(E,ModelInputs,scorefuncs,numTrialPttrn,priorWeightPttrn)
    [~, tempSores, tempModels, allLogEvidence] = wsbmLooper_2(inputData, ...
        modelInputs, ...
        scoreFunc,...
        numTrialPttrn,...
        priorWeightPttrn, ...
        idx, outputDir);
    
    allFitModels(:,idx) = tempModels(:) ;

    allFitLogEvidence(idx, :, :) = allLogEvidence;
    %save the results 
    loopResults(:,idx + 1) = tempSores ;
    
    %TMP
    ppm1.increment();
end
