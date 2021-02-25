function [ cn ] = cliqueNumber( A )
%CLIQUENUMBER Find the size of the maximum clique
% Note: works on undirected graphs

% @input A, an NxN adjacency matrix
% @output cn, a scalar of the clique number (size of maximum clique) on the graph

cn = max(cellfun(@length, maximalCliques( A )));

end

