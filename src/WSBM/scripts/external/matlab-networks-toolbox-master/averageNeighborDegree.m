function av=averageNeighborDegree(A)
%AVERAGENEIGHBORDEGREE Compute the average degree of neighboring nodes for every vertex.
% Note: Works for weighted degrees (graphs) also.
%
% @input A, NxN adjacency matrix
% @output av, a 1xN vector of average neighbor degree
%
% Other routines used: degrees.m, kneighbors.m

% Updated : Changed for function name consistency

% IB: last updated, 3/23/14


av=zeros(1,length(A));   % initialize output vector
[deg,~,~]=degrees(A);

for i=1:length(A)  % across all nodes
  
  neigh=kneighbors(A,i,1);  % neighbors of i, one link away
  if isempty(neigh)
      av(i)=0; 
      continue; 
  end
  av(i)=sum(deg(neigh))/deg(i);
  
end