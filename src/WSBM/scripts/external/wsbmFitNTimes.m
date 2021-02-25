function [ allModels ] = wsbmFitNTimes( adjMat, rStruct , modelInputs , numFits , numCores)
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
    disp('will not run parallel')
    numCores = 0;
end

tempModelStruct = struct() ; 

if numCores > 1
    
    if isempty(gcp)
        parpool(numCores);
    end
    
    %fit numFits times at each iter
    parfor idx=1:numFits % model fits per iteration 

        [~,tempModel] = wsbm(adjMat, ...
            rStruct, ...
            modelInputs{:} ) ;

        % to save space
        tempModel = rmfield(tempModel,'Data') ;
        tempModelStruct(idx).Model = tempModel ;

    end %iterating over number of model fits 
    
else

    %fit numFits times at each iter
    for idx=1:numFits % model fits per iteration 

        [~,tempModel] = wsbm(adjMat, ...
            rStruct, ...
            modelInputs{:} ) ;

        % to save space
        tempModel = rmfield(tempModel,'Data') ;
        tempModelStruct(idx).Model = tempModel ;

    end %iterating over number of model fits

end

allModels = tempModelStruct ;
