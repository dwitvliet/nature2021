function [C1,C2, C] = clustCoeff(A)
%CLUSTCOEFF Compute two clustering coefficients, based on triangle motifs count and local clustering
% C1 = number of triangle loops / number of connected triples
% C2 = the average local clustering, where Ci = (number of triangles connected to i) / (number of triples centered on i)
% Ref: M. E. J. Newman, "The structure and function of complex networks"
% Note: Valid for directed and undirected graphs
%
% @input A, NxN adjacency matrix
% @output C1, a scalar of the average clustering coefficient (definition 1). 
% @output C2, a scalar of the average clustering coefficient (definition 2). 
% @output C, a 1xN vector of clustering coefficients per node (where mean(C) = C2).  
%
% Other routines used: degrees.m, isDirected.m, kneighbors.m, numEdges.m, subgraph.m, loops3.m, numConnTriples.m

% Updated: Returns C vector of clustering coefficients.

% IB, Last updated: 3/23/2014
% Input [in definition of C1] by Dimitris Maniadakis.



n = length(A);
A = A>0;  % no multiple edges
[deg,~,~] = degrees(A);
C=zeros(n,1); % initialize clustering coefficient

% multiplication change in the clust coeff formula
coeff = 2;
if isDirected(A); coeff=1; end

for i=1:n
  
  if deg(i)==1 || deg(i)==0; C(i)=0; continue; end

  neigh=kneighbors(A,i,1);
  edges_s=numEdges(subgraph(A,neigh));
  
  C(i)=coeff*edges_s/(deg(i)*(deg(i)-1));

end

C1=3*loops3(A)/(numConnTriples(A)+2*loops3(A));
C2=sum(C)/n;