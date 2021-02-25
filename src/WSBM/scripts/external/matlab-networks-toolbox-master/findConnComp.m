function comp_mat = findConnComp(A)
%FINDCONNCOMP Find connected components in a graph
% Note: Valid for undirected graphs only
%
% @input A, an NxN adjacency matrix
% @output comp_mat, a cell array where each element is a list of nodes in
% the component.

% Other routines used: findConnCompI.m, degrees.m

% Update: minor optimizations, documentation format
% IB: last updated, 3/23/14

[deg,~,~]=degrees(A);            % degrees
comp_mat={};                       % initialize components matrix

for i=1:length(deg)
    if deg(i)>0
        done=0;
        for j=1:length(comp_mat)
            if any(comp_mat{j}==i)   % i in comp_mat(x).mat
                done=1;
                break
            end
        end
        if not(done)
            comp=findConnCompI(A,i)';
            comp_mat{length(comp_mat)+1}=comp; %#ok<AGROW>
        end
        
    elseif deg(i)==0
        comp_mat{length(comp_mat)+1}=i; %#ok<AGROW>
    end    
end