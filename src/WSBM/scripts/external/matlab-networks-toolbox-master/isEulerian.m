function response=isEulerian(A)
%ISEULERIAN Check if a graph is Eulerian, i.e. it has an Eulerian circuit

% @input A, NxN adjacency matrix
% @output response, boolean response: [0, 1]

% "A connected undirected graph is Eulerian if and only if every graph vertex has an even degree."
% "A connected directed graph is Eulerian if and only if every graph vertex has equal in- and out- degree."
% Note 1: Assume that the graph is connected.
% Note 2: If there is an Eulerian trail, it is reported.
%
% Other routines used: degrees.m, isDirected.m

% Updated: fprintf
% IB: last updated, 3/24/14

response=false;

[degs,indeg,outdeg]=degrees(A);
odd=find(mod(degs,2)==1);

if not(isDirected(A)) & isempty(odd) % if undirected and all degrees are even
  response=true;

elseif isDirected(A) & indeg==outdeg % directed and in-degrees equal out-degrees
  response=true;

elseif numel(odd)==2
  fprintf('isEulerian.m: There is an Eulerian trail from node %2i to node %2i\n',odd(1),odd(2));

end