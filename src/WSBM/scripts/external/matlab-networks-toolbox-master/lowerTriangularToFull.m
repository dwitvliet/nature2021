function [ A_full ] = lowerTriangularToFull( A )
%LOWERTRIANGULARTOFULL Convert a lower triangular matrix to full by copying about the diagonal
%   This method is intended to convert half-adjacency matrices where
%   certain algorithms are easier to implement on lower or upper triangular
%   matrices (avoiding symmetric doubles)

%   Note: method returns 'A' unchanged when A not lower triangular

% @input A, [NxN] adjacency matrix 

% @output A_full, a matrix with lower triangular elements copied across the diagonal
    if(isLowerTriangular(A))
        A_full = tril(A) + tril(A, -1)'; % note, does not copy diagonal twice
    else
        A_full = A;
    end
end

