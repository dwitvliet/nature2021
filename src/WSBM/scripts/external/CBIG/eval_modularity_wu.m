function [ Q , Qvec , Qvec_prcnt ] = eval_modularity_wu(CIJ,ca)
% get the modularity of each commmunity, should sum to traditional Q metric
% https://www.nature.com/nature/journal/v433/n7028/full/nature03288.html
% https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2955452/ eq (8)

% make ca column vec
if ~iscolumn(ca)
   ca = ca'; 
end

% and autofix mat
CIJ = weight_conversion(CIJ,'autofix');

num_coms = length(unique(ca));
tot_sum = sum(sum(CIJ)) ;

% initialize output vector
Qvec = zeros([num_coms 1]) ;

comMat = get_block_mat(CIJ,ca);

for idx=1:num_coms
   
    within_block = comMat(idx,idx) ;
    all_block = sum(comMat(idx,:));
    
    Qvec(idx) = ( within_block / tot_sum ) ...
        - ( all_block / (2 * tot_sum) )^2 ;
    
end

Q = sum(Qvec) ;
Qvec_prcnt = Qvec ./ Q ;

















