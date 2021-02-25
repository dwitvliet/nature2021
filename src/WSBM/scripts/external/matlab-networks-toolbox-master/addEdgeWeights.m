function elc=addEdgeWeights(el)
%ADDEDGEWEIGHTS Add multiple edges in an edge list
% 
% @input el, original (non-compact) edge list
% @output elc, final compact edge list (no row repetitions)
%
% Example: [1 2 2; 2 2 1; 4 5 1] -> [1 2 3; 4 5 1]

% Updated: fixed visited vector for MATLAB.
% IB: last updated, 3/23/2014

el2=[el(:,1), el(:,2)]; % make the edge list searchable w/o the weights
visited=[0, 0];             % mark visited edges

elc=[];
for e=1:size(el,1)
    if sum(ismember(visited,el2(e,:),'rows'))==0  % if not visited yet
        ind=ismember(el2,el2(e,:),'rows');
        ind=find(ind==1);     % these are all the ocurrences of el(e,:)
        elc=[elc; el(e,1), el(e,2), sum(el(ind,3))];
        visited=[visited; el2(e,:)];
    end
end