function [ ret ] = upperTriangularToFull( A )
%UPPERTRIANGULARTOFULL Convert a lower triangular matrix to full by copying about the diagonal
%   This method is intended to convert half-adjacency matrices where
%   certain algorithms are easier to implement on lower or upper triangular
%   matrices (avoiding symmetric doubles)

%   Note: method returns 'A' unchanged when A not upper triangular

% @input A, [NxN] adjacency matrix 

% @output A_full, a matrix with upper triangular elements copied across the diagonal
    if(isUpperTriangular(A))
        ret = triu(A) + triu(A, 1)';
    else
        ret = A;
    end
end

