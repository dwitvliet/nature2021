function [l, path_lengths]  = averagePathLength(A)
%AVERAGEPATHLENGTH Compute average path length for a network - the average shortest path
% Note: works for directed/undirected networks 
%
% @input A, an NxN adjacency matrix
% @output l, a scalar of the average shortest path length
% @output path_lengths, a 1xN vector of shortest path lengths for each node

% INPUTS: adjacency (or weights/distances) matrix, nxn
% OUTPUTS: average path length, optional matrix of pairwise path lengths.
%
% Other routines used: simpleDijkstra.m 

% updated: Handled case where nodes are unreachable t.f. value goes to Inf.

% IB: Last updated, 3/23/14


n=size(A,1);

d = zeros(n);

for i=1:n 
    d(:, i) = simpleDijkstra(A, i);
end
path_lengths = d;
d = d(:);
d = d(d ~= Inf & d ~= 0); % not unreachable and not self
l = nanmean(d); % sum and average across everything but the diagonal