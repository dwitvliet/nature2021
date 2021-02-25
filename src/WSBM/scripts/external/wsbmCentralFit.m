function [ centralModel , allModels ] = wsbmCentralFit( adjMat, rStruct , modelInputs , numFits , priorMu, numTrialPttrn, priorWeightPttrn)
% will fit wsbm specificed number iterations
% and will return 'most central' based on variation of information between
% community labels
% 
% user has the option to also provide prior labels to additionally multiply
% by when computing the central mode. this will be useful if you want to
% keep the central model 'near' a prior 

% edit... numFits can be a vector, specifying num fits at each level of vec

if nargin < 4
    disp('need more args')
    return
end

if nargin < 5
    disp('no prior read, this is cool')
    priorMu = '';
end

if ~exist('numTrialPttrn','var') || isempty(numTrialPttrn)
    numTrialPttrn = [] ;
end

if ~exist('priorWeightPttrn','var') || isempty(priorWeightPttrn)
    priorWeightPttrn = [] ;
end

tempModelStruct = struct() ; 

%fit numFits times at each iter
for idx=1:numFits % model fits per iteration 
    
    if isempty(numTrialPttrn)
        
        [~,tempModel] = wsbm(adjMat, ...
            rStruct, ...
            modelInputs{:} ) ;
    else
        
        if numel(rStruct) == 1
            r_size = rStruct ;
        else
            r_size = size(rStruct,1) ;
        end
             
       %begin with uninformative prior
       mu_prior = make_WSBM_prior(...
           ones(r_size,size(adjMat,1)), 0) ;
        
       for jdx=1:numel(numTrialPttrn)
        
            % the parameters set after 'modelInputs' will take precendent
            % over anything specified in that cell array
            [~,tempModel] = wsbm(adjMat, ...
                rStruct, ...
                modelInputs{:},...
                'numTrials', numTrialPttrn(jdx), ...
                'mu_0', mu_prior ) ;
                         
            if isempty(priorWeightPttrn)
                mu_prior = make_WSBM_prior(Model,jdx) ;
            else
                mu_prior = make_WSBM_prior(Model,priorWeightPttrn(jdx));
            end 
       end
           
    end

    % to save space
    tempModel = rmfield(tempModel,'Data') ;
    tempModelStruct(idx).Model = tempModel ;

end %iterating over number of model fits

centralModel = central_model( tempModelStruct , priorMu) ;
allModels = tempModelStruct ;
