function d = simpleDijkstra(A,s)
%SIMPLEDIJKSTRA Dijkstra's shortest path algorithm
% Returns the distances from a single vertex to all others, doesn't save the path
%
% @input A, an NxN adjacency matrix
% @input s, a scalar of the start node
% @output d, a 1xn vector of shortest path lengths from 's' to all other
% nodes.

% Note: Works for a weighted/directed graph.

% Updated 3/8/14: increased efficency 10x by unrolling for loop, removing
% setdiff.

% IB: Last updated, 3/23/14

n=length(A);
d = inf*ones(1,n); % distance s-all nodes
d(s) = 0;    % s-s distance
T = 1:n;    % node set with shortest paths not found yet

while not(isempty(T))
    [~,ind] = min(d(T));
    
    idx_update = find(A(T(ind), T) > 0 & d(T) > d(T(ind)) + A(T(ind), T));
    if(~isempty(idx_update))
        d(T(idx_update)) = d(T(ind)) + A(T(ind), T(idx_update));
    end    
    T(ind) = [];
end