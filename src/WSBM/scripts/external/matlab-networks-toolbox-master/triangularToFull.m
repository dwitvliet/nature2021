function [ ret ] = triangularToFull( A )
%TRIANGULARTOFULL Convert a triangular matrix (lower or upper) to full by copying about the diagonal
%   This method is intended to convert half-adjacency matrices where
%   certain algorithms are easier to implement on lower or upper triangular
%   matrices (avoiding symmetric doubles)

%   Note: method returns 'A' unchanged when A not triangular (lower or upper)

% @input A, [NxN] adjacency matrix 

% @output A_full, a matrix with lower or upper triangular elements copied across the diagonal


if(isLowerTriangular(A))
    ret = lowerTriangularToFull(A);
elseif(isUpperTriangular(A))
    ret = upperTriangularToFull(A);
else
    ret = A;
end

