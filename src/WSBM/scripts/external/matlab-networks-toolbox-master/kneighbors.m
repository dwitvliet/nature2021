function neighbors = kneighbors(A,ind,k)
%KNEIGHBORS Finds the adjacency list of nodes at distance k from 'ind'
%
% @input A, NxN adjacency matrix
% @input ind, a scalar node # of the search index
% @input k, a scalar for distance to search
% @output neighbors, a adjacency list of nodes reachable from 'ind' in 'k' hops 
% INPUTS: adjacency matrix (nxn), start node index, k - number of links
% OUTPUTS: vector of k-neighbors indices
%
% Updated: For readability.

% IB: last updated, 3/23/14

adjk = A;
for i=1:k-1 
    adjk = adjk*A; 
end;

neighbors = find(adjk(ind,:)>0);