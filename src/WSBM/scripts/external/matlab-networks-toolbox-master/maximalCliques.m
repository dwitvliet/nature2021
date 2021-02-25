function [ MC ] = maximalCliques( A )
%MAXIMALCLIQUES Find maximal cliques using the Bron-Kerbosch algorithm
%   Given a graph's boolean adjacency matrix, A, find all maximal cliques 
%   on A using the Bron-Kerbosch algorithm in a recursive manner.  The 
%   graph is required to be undirected and must contain no self-edges.

%   This script is adapted from Jeffery Wildman's script
%   (http://www.mathworks.com/matlabcentral/fileexchange/30413-bron-kerbosch-maximal-clique-finding-algorithm),
%   implementing Adrian Wilder's suggested speed-up (a factor of > 200 speed up on my test data)

%   **For speed (e.g. I am using this as a subroutine on >1M small subgraphs), there is no input checking, use standard adjacency matrices (undirected, no loops) **

%   @input A, NxN boolean adjacency matrix
%   @output MC, K-length cell array of nodes adjacent in each of K maximal cliques on the graph.


%   Ref: Bron, Coen and Kerbosch, Joep, "Algorithm 457: finding all cliques
%   of an undirected graph", Communications of the ACM, vol. 16, no. 9, 
%   pp: 575–577, September 1973.
%
%   Ref: Cazals, F. and Karande, C., "A note on the problem of reporting 
%   maximal cliques", Theoretical Computer Science (Elsevier), vol. 407,
%   no. 1-3, pp: 564-568, November 2008.

%   IB: last updated, 3/23/14

%% preallocate
n = size(A,2); % number of vertices
MC = cell(n,1); % storage for maximal cliques
R = false(n,1); % currently growing clique
P = true(n,1); % prospective nodes connected to all nodes in R
X = false(n,1); % nodes already processed
iclique=0;
A=A.'; %this speeds up some of the calculations below because we do not have to transpose A for each recursion

%% run
bron_kerbosch(R,P,X);

%% trim output
MC((iclique+1):end)=[];

    function [] = bron_kerbosch( R, P, X )
        
        if ~any(P | X)
            % report R as a maximal clique
            iclique=iclique+1;
            MC{iclique}=find(R);
        else
            % choose pivot
            ppivots = P | X; % potential pivots
            binP = zeros(1,n);
            binP(P) = 1; % binP contains ones at indices equal to the values in P
            pcounts = binP*double(A(:,ppivots)); % cardinalities of the sets of neighbors of each ppivots intersected with P
            % select one of the ppivots with the largest count
            [~,ind] = max(pcounts);
            temp_u=find(ppivots,ind,'first');
            u_p=temp_u(ind);
            
            for u = find(~A(:,u_p) & P).' % all prospective nodes who are not neighbors of the pivot
                P(u)=false;
                Rnew = R;
                Rnew(u)=true;
                Nu = A(:,u);
                Pnew = P & Nu;
                Xnew = X & Nu;
                bron_kerbosch(Rnew, Pnew, Xnew);
                X(u)=true;
            end
        end
        
    end % BKv2
end % maximalCliques
