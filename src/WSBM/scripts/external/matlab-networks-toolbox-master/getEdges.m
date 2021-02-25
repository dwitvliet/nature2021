function edges = getEdges(graph,type)
%GETEDGES Returns the list of edges for graph varying representation types
% @input graph, (multiple formats: adjacency matrix, adjacencylist, edgelist, inc) (see below)
% @input type, a string specifier of the format type
% @output edges, Mx3 edge list where third column is weight.
% 
% Note 1: 'type' can be: 'adj','edgelist','adjlist' (neighbor list), 'inc' (incidence matrix)
% Note 2: symmetric edges will appear twice, also in undirected graphs, (i.e. [n1,n2] and [n2,n1])
% Other routines used: adj2edgeL.m, adjL2edgeL.m, inc2edgeL.m
%
% Example representations of a directed triangle: 1->2->3->1
%           'adj' - [0 1 0; 0 0 1; 1 0 0]
%           'adjlist' - {1: [2], 2: [3], 3: [1]}
%           'edgelist' - [1 2; 2 3; 3 1] or [1 2 1; 2 3 1; 3 1 1] (1 is the edge weight)
%           'inc' - [-1  0  1
%                     1 -1  0
%                     0  1 -1]
%
% Updated:  fprintf
% IB: last updated, 3/24/14



if strcmp(type,'adj')
    edges=sortrows(adj2edgeL(graph));
    
elseif strcmp(type,'edgelist')
    edges=graph; % the graph structure is the edge list
    
elseif strcmp(type,'adjlist')
    edges=sortrows(adjL2edgeL(graph));
    
elseif strcmp(type,'inc')
    edges=sortrows(inc2edgeL(graph));
else
    fprintf('ERROR: "type" input can only be "adj" (adjacency, nxn matrix), "edgelist" (mx3 matrix)\n, "adjlist" (neighbor list, nx1 cell) and "inc" incidence (nxm matrix)\n')
end