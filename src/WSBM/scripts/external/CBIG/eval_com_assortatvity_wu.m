function [ comMaxAssorVec , nodMaxAssorVec ] = eval_com_assortatvity_wu(CIJ,ca) 
% return the maximum assortativity from Betzel
% https://arxiv.org/pdf/1702.02807.pdf

% make ca column vec
if ~iscolumn(ca)
   ca = ca'; 
end

% and autofix mat
CIJ = weight_conversion(CIJ,'autofix');

% number of communities
num_coms = length(unique(ca));

% number of nodes
num_nodes = size(CIJ,1);

% initialize community, size mat
comMat = zeros(num_coms) ;
sizeMat = zeros(num_coms) ;

% initialize output vecs
comMaxAssorVec = zeros([ num_coms 1]);
nodMaxAssorVec = zeros([ num_nodes 1]);

% community indicies
caIdx = ~~dummyvar(ca);

for idx=1:num_coms
    for jdx=1:num_coms
       
        comMat(idx,jdx) = sum(sum(CIJ(caIdx(:,idx),caIdx(:,jdx))));
    
        % to match paper, calc sizes
        sizeMat(idx,jdx) = numel(CIJ(caIdx(:,idx),caIdx(:,jdx)));
        
    end
end

% find the max assort at each community 
% eq (7) https://arxiv.org/pdf/1702.02807.pdf
for idx=1:num_coms
    
    ind = ones([num_coms 1]);
    ind(idx) = 0 ;
    betwn = comMat(idx,~~ind) ./ sizeMat(idx,~~ind);
    
    maxBetwn = max(betwn);   
    within = comMat(idx,idx) ./ sizeMat(idx,idx);
    
    comMaxAssorVec(idx) = within - maxBetwn ;
    
end

% find assortatvitiy at each node
% eq (8) https://arxiv.org/pdf/1702.02807.pdf
for idx=1:num_nodes
      
    nodComMat = zeros([num_coms 1]) ;
    
    % measure this nodes weighted density to each community
    % normalized by size of that community
    for jdx=1:num_coms
        
        nodComMat(jdx) = sum(CIJ(idx,ca == jdx)) / sum(ca == jdx); %sizeMat(jdx,jdx);
        
    end
    
    nodeCom = ca(idx);
    nodWithin = nodComMat(nodeCom) ;
    
    indBetwn = ones([num_coms 1]);
    indBetwn(nodeCom) = 0 ;
    
    nodeMaxBetwn = max(nodComMat(~~indBetwn));
    
    nodMaxAssorVec(idx) = nodWithin - nodeMaxBetwn ;
    
end










