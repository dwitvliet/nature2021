function [ num_cc, max_cc ] = connectedComponentMeasures( A )
%CONNECTEDCOMPONENTMEASURES find the size of the largest connected component
% This method finds the size of the largest connected component, a
% frequently used network measurement/feature.

% @input A, Adjacency matrix
% @output num_cc a scalar value of the number of connected components
% @output max_cc a scalar value of the size (# of nodes) of the largest connected component

% IB: last updated, 3/23/14

cc = cellfun(@length, findConnComp(A));
max_cc = max(cc);
num_cc = length(cc);


end

