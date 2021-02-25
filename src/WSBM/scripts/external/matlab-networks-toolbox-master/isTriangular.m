function [ ret ] = isTriangular( A, enforce_strictly )
%ISTRIANGULAR Tests if a matrix is triangular (upper or lower)
%   This method is intended to check half-adjacency matrices where
%   certain algorithms are easier to implement on lower or upper triangular
%   matrices (avoiding symmetric doubles)

% @input A, [NxN] adjacency matrix 
% @input enforce_strictly[optional, default: false], a boolean choice to
% ensure no nonzero elements are along the diagonal. 

% @output boolean response if A is triangular

if(~exist('enforce_strictly', 'var') || isempty(enforce_strictly) || enforce_strictly == 0)
    enforce_strictly = false(1);
else
    enforce_strictly = true(1);
end

ret = isUpperTriangular(A, enforce_strictly) | isLowerTriangular(A, enforce_strictly);
end

