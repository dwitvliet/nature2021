function [ ret ] = isLowerTriangular( A, enforce_strictly )
%ISLOWERTRIANGULAR Tests if a matrix is lower-triangular
%   This method is intended to check half-adjacency matrices where
%   certain algorithms are easier to implement on lower or upper triangular
%   matrices (avoiding symmetric doubles)

% @input A, [NxN] adjacency matrix 
% @input enforce_strictly[optional, default: false], a boolean choice to
% ensure no nonzero elements are along the diagonal. 

% @output boolean response if A is lower triangular

if(~exist('enforce_strictly', 'var') || isempty(enforce_strictly) || enforce_strictly == 0)
    enforce_strictly = false(1);
else
    enforce_strictly = true(1);
end

a = any(any(triu(A, 1)));
b = any(diag(A));
if(enforce_strictly && ~b && ~a)
    ret = 1;
elseif(~enforce_strictly && ~a)
    ret =1;
else
    ret = 0;
end

