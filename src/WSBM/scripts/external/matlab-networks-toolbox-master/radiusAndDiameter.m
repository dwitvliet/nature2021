function [r, d] = radiusAndDiameter(A)
%RADIUSANDDIAMETER calculate the radius and diameter of the graph

% This method returns radius the minimum graph eccentricity (Minimum
% longest shortest path), and the maximum graph accentricity (maximum longest shortest path).
%
% @input adj, NxN adjacency matrix
% @output r, a scalar reporting the minimum longest shorest path
% @output d, a scalar reporting the maximum longest shorest path
%
% Other routines used: simpleDijkstra.m
% IB: last updated, 3/23/14

rs = zeros(size(A, 1), 1);
for i=1:size(A,1)
    d=simpleDijkstra(A,i);
    idx_good = d ~= Inf & d ~= 0;
    if(any(idx_good))
        rs(i) = max(d(idx_good));
    else
        rs(i) = NaN;
    end
end
r = nanmin(rs);
d = nanmax(rs);