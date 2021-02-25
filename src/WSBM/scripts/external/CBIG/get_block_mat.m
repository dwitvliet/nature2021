function [weiBM,avgWeiBM,binBM,avgBinBM,stdWeiBM,sizesMat] = get_block_mat(CIJ,ca,excludeNaN)
% given an adjacency matrix + community affiliations, return a block matrix
%
% inputs: 
%           CIJ:        adjacency matrix
%           ca:         community affiliation vector
%
%           excludeNaN: if you dont exlude NaN, when calculating the
%           averages we will count nan edges still in the average... else,
%           if excludeNaN == 1, treat the NaNs as throwaway
%
% returns: 
%           weiBM:      sum of weights block matrix
%           avgWeiBM:   average weighted block matrix
%           binBM:      sum of binary weights block matrix
%           avgBinBM:   average binary block matrix
%           stdWeiBM:   std weights block matrix 
%           sizesMat:   number of edges per block (incl. NaN edges)
%
% Josh Faskowtiz IU

if nargin < 3
    excludeNaN = 0 ;
end

% make ca column vec
if ~iscolumn(ca)
   ca = ca'; 
end

% number coms
nBlocks = length(unique(ca));

% do this in case the user has provided ca with numeric gaps, eg ca=[1,2,4] 
orderedBlocks = sort(unique(ca));

% number nodes per block
blockSizes = histcounts(sort(ca));

sizesMat = bsxfun(@times,...
    (bsxfun(@times,ones(nBlocks),blockSizes)),...
    blockSizes');

% initialize outputs
weiBM = zeros([nBlocks nBlocks]);
avgWeiBM = zeros([nBlocks nBlocks]);
binBM = zeros([nBlocks nBlocks]);
avgBinBM = zeros([nBlocks nBlocks]);
stdWeiBM = zeros([nBlocks nBlocks]);

% now loop it
for idx = 1:nBlocks % rows
    for jdx = 1:nBlocks %columns

        tmp = CIJ(ca == orderedBlocks(idx),ca == orderedBlocks(jdx));
        
        % weighted
        weiBM(idx,jdx) = nansum(tmp(:));
        avgWeiBM(idx,jdx) = nanmean(tmp(:));
        
        % binary
        tmpbin = (isnan(tmp) .* tmp) + (tmp > 0) ;      
        binBM(idx,jdx) = nansum(tmpbin(:));
        avgBinBM(idx,jdx) = nanmean(tmpbin(:));
        
        % std
        stdWeiBM(idx,jdx) = nanstd(tmp(:));
         
    end
end

% if we are not excluding the NaN's when averaging, we can just divide the
% avg block mats by the size Mat, and overwrite what was previously
% computed
if excludeNaN == 0
    avgWeiBM = weiBM ./ sizesMat;
    avgBinBM = binBM ./ sizesMat;
end


