function [ templateMatBothHemi , meanLens , meanCoorMM, meanCountMat ] = make_template_mat(dataStruct, leftNodes, rightNodes, initThresh)
%% make average mat, and return it
% we will use function that preserves mat lengths
% therefore this function needs the additional function: fcn_group_average

selectNodes = [ leftNodes rightNodes ] ;

%% make average distCoorMM
% will compute both for fun but will only use the lengths of the fibers
% for when we run it, but we can change this

% euclidean distance

arrayCoorMM = zeros([size(dataStruct(1).distCoorMM) length(dataStruct)]) ;
for i=1:length(dataStruct)
    arrayCoorMM(:,:,i) = dataStruct(i).distCoorMM ; 
end
meanCoorMM = mean(arrayCoorMM,3) ;
% 
% clear arrayCoorMM

% tract distance on average

arrayLens = zeros([size(dataStruct(1).lensMat) length(dataStruct)]) ;
for i=1:length(dataStruct)
    arrayLens(:,:,i) = dataStruct(i).lensMat ;
end
meanLens = mean(arrayLens,3) ;

%% make average mat
% by looping over data given 

% pre-allocate the matrix to store subject-level mats
sizeAllNodes = length(selectNodes) ;
arrayWeightMats = zeros(sizeAllNodes, sizeAllNodes, length(dataStruct)) ;
arrayRawCountMats = zeros(sizeAllNodes, sizeAllNodes, length(dataStruct)) ;

for idx=1:length(dataStruct)
           
    tmp = dataStruct(idx).countVolNormMat(selectNodes,selectNodes) ;
    
    %get size of square mat
    n=size(tmp,1);
    %make zero across diagonal
    tmp(1:n+1:end) = 0;
    
    %make a mask based on thresh
    temp_mask = dataStruct(idx).countMat(selectNodes,selectNodes) > initThresh ;    
    temp_mask(temp_mask > 0) = 1 ;
    
    % threshold temp matrix by mask
    tmp = tmp .* temp_mask ;
    arrayWeightMats(:,:,idx) = tmp ;
    
    arrayRawCountMats(:,:,idx) = dataStruct(idx).countMat(selectNodes,selectNodes) ...
        .* temp_mask;
    
end

% make the lh rh membership mat
hemi_id = ones(sizeAllNodes,1);
tmp = length(leftNodes) + 1 ; 
hemi_id(tmp:end) = 2 ;

% using mean lengths here
template_mask = fcn_group_average(arrayWeightMats, ...
    meanLens(selectNodes,selectNodes), ...
    hemi_id) ; 

% get number of edge existences
numEdg = sum(arrayWeightMats>0, 3);

avrWei = sum(arrayWeightMats, 3)./ numEdg;
avrWei(numEdg == 0) = 0;  % remove NANs

%mask the avrWei mat by the threshold mat
templateMatBothHemi = avrWei .* template_mask ;

%also return count mat 
meanCountMat = (sum(arrayRawCountMats,3) ./numEdg) .* template_mask ;
meanCountMat(numEdg == 0) = 0;

