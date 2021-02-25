function [ Model ] = wsbmFitWPttrn( adjMat, rStruct , modelInputs , initMu, numTrialPttrn, priorWeightPttrn)
% will fit wsbm specificed number iterations pattern and prior weight
% pattern
% if you want the wsbm fit with no prior, just set initMu to a uniform
% prior and it will have no effects 

if nargin < 5
    disp('need more args')
    return
end

if ~exist('priorWeightPttrn','var') || isempty(priorWeightPttrn)
    priorWeightPttrn = [] ;
end

mu_prior = initMu ;

for idx=1:length(numTrialPttrn)

    % the parameters set after 'modelInputs' will take precendent
    % over anything specified in that cell array
    [~,tmpModel] = wsbm(adjMat, ...
        rStruct,...
        modelInputs{:},...
        'numTrials', numTrialPttrn(idx), ...
        'mu_0', mu_prior ) ;

    if isempty(priorWeightPttrn)
        mu_prior = make_WSBM_prior(tmpModel,idx) ;
    else
        mu_prior = make_WSBM_prior(tmpModel,priorWeightPttrn(idx));
    end

end

Model = tmpModel ;
