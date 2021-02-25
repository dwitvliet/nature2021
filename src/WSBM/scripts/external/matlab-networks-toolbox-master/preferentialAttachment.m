function el = preferentialAttachment(n,m)
%PREFERENTIALATTACHMENT simple preferential attachment for network growth.

% The probability that a new vertex attaches to a given old vertex
%                           is proportional to the (total) vertex degree.

% @input n final number of desired vertices
% @input m # of edges to attach at each step
% @output el, edge list in Mx3 format. 

% Note 1: Vertices arrive one at a time.
% Note 2: Assume undirected simple graph.
% Source: Newman, "The Structure and Function of Complex Networks"
%         B-A., "Emergence of Scaling in Random Networks"

% Other routines used: weightedRandomSample.m
% Updated: fprintf, documentation

% IB: last updated, 3/24/14
%##################################################################



vertices = 2;
if not(vertices<=n); fprintf('Specify more than 2 nodes.\n');  return; end
el=[1 2 1; 2 1 1];      % start with one edge


while vertices < n
  
  vertices=vertices+1;  % add new vertex

  if m>=vertices
    for node=1:vertices-1
      el = [el; node vertices 1];
      el = [el; vertices node 1];  % add symmetric edge
    end
    continue
  end
    
  deg=[];        % compute nodal degrees for this iteration
  for v=1:vertices; deg=[deg; numel(find(el(:,1)==v))]; end
  

  % add m edges
  r = weightedRandomSample(m,[1:vertices],deg/sum(deg));
  while not(length(unique(r))==length(r))
    r = weightedRandomSample(m,[1:vertices],deg/sum(deg));
  end
  
  for node=1:length(r)
    el = [el; r(node) vertices 1];
    el = [el; vertices r(node) 1];      
  end      
  
end