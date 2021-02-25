function [ d ] = edgeDensity( A )
%EDGEDENSITY An alias method for linkDensity
% Edge density is defined as the number of edges divided by 
% number_of_nodes(number_of_nodes-1)/2 where the latter is the maximum possible number of edges.
%
% @input A, NxN adjacency matrix
% @output d, scalar density value between 0 and 1.
%
% Note 1: The graph has to be non-trivial (more than 1 node).
% Note 2: Routine works for both directed and undirected graphs.
%
% Other routines used: numNodes.m, numEdges.m, isDirected.m

% Updated: documentation, added edgeDensity alias.  
% IB: last update 3/24/14

n = numNodes(A);

coeff = 2;
if isDirected(A); coeff = 1; end

d = coeff*numEdges(A)/(n*(n-1));
end

