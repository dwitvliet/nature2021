function [ centralModel , simVec , simTriu ] = central_model( modelsStruct , priorMu , priorWeightTune)
% function for iterating over a stuct to find the central model, in terms
% of community labels, mu.  Have the option of providing a prior too
% edited this to also take in cell array of models

if nargin < 2
    disp('no prior read, this is cool')
    priorMu = '';
end

if ~exist('priorWeightTune','var') || isempty(priorWeightTune)
    priorWeightTune = -1;
end

if iscell(modelsStruct)
    modelsStruct = cell2struct(modelsStruct(:),{'Model'},2) ;
end

% preallocate 
numModels = length(modelsStruct) ;
simMat = zeros([ numModels numModels ]) ;

for idx=1:numModels
    for jdx=1:numModels
    
        % for each model in the stuct, get harsh mu, which we will take to be the 
        % finte communty strucutre
        [ ~ , mu1 ] = make_WSBM_prior(modelsStruct(idx).Model) ;
        [ ~ , mu2 ] = make_WSBM_prior(modelsStruct(jdx).Model) ;

        % use negative varation of information so we can maximize
        % it later, and then sum across one dimension 
%         % (the matrix we summing over will be symetric
% disp(mu1);
% disp('ha');
% disp(mu2);
        simMat(idx,jdx) = -varInfo( mu1, mu2) ;

    end
end

triuIdx = triu(ones(numModels)) ;
simTriu = simMat(~~triuIdx);
simVec = sum(simMat,2) ;

% check to see if all resulting mu_0 are the same....
totSim = sum(simVec) ; 
if totSim == 0
    % if all the fits are the same, just return the first as the
    % central model
    centralModel = modelsStruct(1).Model;
    % return from function 
    return
else
    %this should be then normal case
    
    %check if there was a priorMu
    % if this priorMu is not empty 
    if ~isempty(priorMu)
        
        disp('adding prior dist as well')
        
        viToPrior = zeros([ numModels 1 ]);
        
        for idx=1:numModels 

            % hash prior again
            [ ~ , mu ] = make_WSBM_prior(modelsStruct(idx).Model) ;
            viToPrior(idx,1) = -varInfo(mu,priorMu) ;

        end
        
        if priorWeightTune < 0
            [ ~ , idxStruct ] = max(sum( (simVec + viToPrior), 2)) ; 
        else
            [ ~ , idxStruct ] = max(sum( simVec .* (1-priorWeightTune) ...
                + (viToPrior .* priorWeightTune) ) , 2) ;
        end
        
    else
        % no prior, just return the central model
        [ ~ , idxStruct ] = max(sum(simVec, 2)) ; 
        
    end
    
    %return model
    centralModel = modelsStruct(idxStruct).Model ;

end
