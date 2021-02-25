function A_sp = minSpanTree(A)
%MINSPANTREE Prim's minimal spanning tree algorithm

% @input A, NxN adjacency matrix
% @output A_sp, NxN adjacency matrix of minimum spanning tree

% Prim's alg idea:
%  start at any node, find closest neighbor and mark edges
%  for all remaining nodes, find closest to previous cluster, mark edge
%  continue until no nodes remain
%
% INPUTS: graph defined by adjacency matrix, nxn
% OUTPUTS: matrix specifying minimum spanning tree (subgraph), nxn
%
% Other routines used: isConnected.m
% Updated: cleaned fprintf

% IB Last updated: 3/24/14


A_sp = [];

% check if graph is connected:
if not(isConnected(A)); fprintf('This graph is not connected. No spanning tree exists.\n'); return; end

n = length(A); % number of nodes
A_sp = zeros(n);   % initialize tree

A(A==0)=inf; % set all zeros in the matrix to inf

conn_nodes = 1;        % nodes part of the min-span-tree
rem_nodes = [2:n];     % remaining nodes

while ~isempty(rem_nodes)
  
  [minlink]=min(min(A(conn_nodes,rem_nodes)));
  ind=find(A(conn_nodes,rem_nodes)==minlink);

  [ind_i,ind_j] = ind2sub([length(conn_nodes),length(rem_nodes)],ind(1));

  i=conn_nodes(ind_i); j=rem_nodes(ind_j); % gets back to adj indices
  A_sp(i,j)=1; A_sp(j,i)=1;
  conn_nodes = [conn_nodes j]; %#ok<AGROW>
  rem_nodes = setdiff(rem_nodes,j);
  
end