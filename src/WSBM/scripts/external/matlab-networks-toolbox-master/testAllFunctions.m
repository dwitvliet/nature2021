% Test code for "Octave tools for Network Analysis"

% IB: 3/24/14, TODO: newmanEigenvectorMethod test fails. Investigate.

% ** start ** 

clear all
close all

% Set of test graphs, in various formats =========
one_directed_edge = [0 1; 0 0];
one_double_edge = [0 2; 2 0]; 
bowtie=[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 1 0 0; 0 0 1 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
disconnected_bowtie =[0 1 1 0 0 0; 1 0 1 0 0 0; 1 1 0 0 0 0; 0 0 0 0 1 1; 0 0 0 1 0 1; 0 0 0 1 1 0];
bowtie_edgeL = [1,2,1; 1,3,1; 2,3,1; 3,4,1; 4,5,1; 4,6,1; 5,6,1];
bowtie_edgeL = sortrows(symmetrizeEdgeL(bowtie_edgeL));
bowtie_edgeL_loop = [bowtie_edgeL; 4 4 1];
bowtie_adjL = {[2,3],[1,3],[1,2,4],[3,5,6],[4,6],[4,5]};
undirected_cherry = [1,2,1; 2,1,1; 1,3,1; 3,1,1];
directed_cherry = [1,2,1; 1,3,1];
undirected_triangle=[0 1 1; 1 0 1; 1 1 0];
undirected_triangle_selfloops = [1 1 1; 1 1 1; 1 1 0];
undirected_triangle_incidence = [1 1 0; 1 0 1; 0 1 1];
directed_triangle=[0 1 0; 0 0 1; 1 0 0];
square = [0 1 0 1; 1 0 1 0; 0 1 0 1; 1 0 1 0];
star = edgeL2adj(canonicalNets(5,'star'));      % adjacency matrix
% ================================================


% Testing getNodes.m =============================
fprintf('testing getNodes.m\n')

assert_mat(getNodes(bowtie,'adj'), [1:6])
N = randi(100);
assert_mat(getNodes(randomDirectedGraph(N),'adj'),[1:N])
assert_mat(getNodes(randomGraph(10),'adj'),[1:10])
assert_mat(getNodes(bowtie_adjL,'adjlist'),[1:6])
assert_mat(getNodes(directed_cherry,'edgelist'),[1:3])
assert_mat(getNodes(undirected_cherry,'edgelist'),[1:3])
assert_mat(getNodes(undirected_triangle_incidence,'inc'),[1:3])
% ================================================


% Testing getEdges.m =============================
fprintf('testing getEdges.m\n')

assert_mat(getEdges(bowtie,'adj'),bowtie_edgeL)
assert_mat(getEdges(bowtie_adjL,'adjlist'),bowtie_edgeL)
assert_mat(getEdges(directed_cherry,'edgelist'),directed_cherry)
assert_mat(getEdges(undirected_cherry,'edgelist'),undirected_cherry)
assert_mat(getEdges(undirected_triangle_incidence,'inc'),[1,2,1; 1,3,1; 2,1,1; 2,3,1; 3,1,1; 3,2,1])
% ================================================

% testing numNodes.m =============================
fprintf('testing numNodes.m\n')

randint = randi(101);
assert_mat(numNodes(randomGraph(randint)),randint)
assert_mat(numEdges(edgeL2adj(directed_cherry)),2)
assert_mat(numNodes(bowtie),6)
% ================================================

% testing numEdges.m =============================
fprintf('testing numEdges.m\n')

assert_mat(numEdges(bowtie),7)
assert_mat( numEdges(undirected_triangle_selfloops), 5 )
assert_mat(numEdges(one_double_edge),2)
assert_mat(numEdges(edgeL2adj(bowtie_edgeL_loop)),8)
% ================================================

% testing linkDensity.m ==========================
fprintf('testing linkDensity.m\n')

randint = randi(101);
assert_mat(linkDensity(edgeL2adj(canonicalNets(randint,'tree',2))),2/randint)
assert_mat(linkDensity(bowtie),2.0*7/(6*5))
% ================================================


% testing selfLoops.m ============================
fprintf('testing selfLoops.m\n')

assert_mat(selfLoops(undirected_triangle_selfloops),2)
assert_mat(selfLoops(bowtie),0)
% ================================================


% testing multiEdges.m ===========================
fprintf('testing multiEdges.m\n')

assert_mat(multiEdges(one_double_edge),2)
assert_mat(multiEdges([0 2 1; 2 0 1; 1 1 0]))  % triangle with one double edge
assert_mat(multiEdges([0 0 1; 2 0 0; 0 1 0]))  % directed triangle with 1 double edge
assert_mat(multiEdges(randomGraph(randi(15))),0)
% ================================================

% testing averageDegree.m ========================
fprintf('testing averageDegree.m\n')

assert_mat(averageDegree(square),2)
assert_mat(averageDegree(bowtie),2+1.0/3)
% ================================================

% testing numConnComp.m ==========================
fprintf('testing numConnComp.m\n')
nc=numConnComp(disconnected_bowtie);
assert_mat(numConnComp(disconnected_bowtie),2)

randint = randi(51);
Adj=zeros(randint*30);
for x=1:randint
  adj=randomGraph(30,0.5);
  Adj(30*(x-1)+1:30*x,30*(x-1)+1:30*x)=adj;
end
assert_mat(numConnComp(Adj),randint)
% ================================================


% testing findConnComp.m =========================
fprintf('testing findConnComp.m\n')

assert_mat(findConnCompI(disconnected_bowtie,1),[1,2,3])
assert_mat(findConnComp(disconnected_bowtie),{[1,2,3],[4,5,6]})

clear modules
modules{1}=[0];
randint = randi(21);
Adj = []; adj = [];

% make up a matrix (Adj) of randint disconnected components (adj)
for x=1:randint
  randsecint = randi(25)+5;
  lastnode = modules{length(modules)}(length(modules{length(modules)}));
  modules{length(modules)+1} = [lastnode+1:lastnode+randsecint]; 
  while isempty(adj) | not(isConnected(adj)) | not(length(adj)==randsecint); adj=randomGraph(randsecint,0.5); end

  Adj(length(Adj)+1:length(Adj)+randsecint,length(Adj)+1:length(Adj)+randsecint)=adj; 
end

modules=modules(2:length(modules));
assert_mat(findConnComp(Adj),modules)
% ================================================


% testing giantComponent.m =======================
fprintf('testing giantComponent.m\n')

clear modules
modules{1}=[0];
randint = randi(20)+1;
Adj = []; adj = [];

% make up a matrix (Adj) of randint disconnected components (adj)
for x=1:randint
  randsecint = randi(25)+5;
  lastnode = modules{length(modules)}(length(modules{length(modules)}));
  modules{length(modules)+1} = [lastnode+1:lastnode+randsecint]; 
  while isempty(adj) | not(isConnected(adj)) | not(length(adj)==randsecint); adj=randomGraph(randsecint,0.5); end
  Adj(length(Adj)+1:length(Adj)+randsecint,length(Adj)+1:length(Adj)+randsecint)=adj; 
end
modules=modules(2:length(modules));
L = [];
for m=1:length(modules); L = [L, length(modules{m})]; end;
[maxL,maxind] = max(L);
assert_mat(giantComponent(Adj), subgraph(Adj,modules{maxind}))
% ================================================


% ================================================
% Testing tarjan.m ===============================
fprintf('testing tarjan.m\n')

L = {}; L{1} = 2; L{2} = 1;
GSCC = tarjan(L);
assert_mat(length(GSCC),1)
assert_mat(GSCC{1},[1,2])

L = {}; L{1} = 2; L{2} = [];
GSCC = tarjan(L);
assert_mat(length(GSCC),2)
assert_mat(GSCC{1},[2])
assert_mat(GSCC{2},[1])


L={}; L{1}=[2,3]; L{2}=[1]; L{3}=[1]; L{4}=[1]; % cherry tree (binary) + extra node
GSCC = tarjan(L);
assert_mat(length(GSCC),2)
assert_mat(GSCC{1},[1,2,3])
assert_mat(GSCC{2},4)


L={}; L{1}=[2,3]; L{2}=[1,3]; L{3}=[1,2]; L{4}=[1]; % triangle with extra node
GSCC = tarjan(L);
assert_mat(length(GSCC),2)
assert_mat(GSCC{1},[1,2,3])
assert_mat(GSCC{2},4)


L={}; L{1}=[2,3]; L{2}=[1,3]; L{3}=[1,2,4]; L{4}=[5,6]; L{5}=[4,6]; L{6}=[4,5];
GSCC = tarjan(L);
assert_mat(length(GSCC),2)
assert_mat(length(GSCC{1}),3)
assert_mat(length(GSCC{2}),3)

L={}; L{1}=[2,3]; L{2}=[1,3]; L{3}=[1,2]; L{4}=[5,6]; L{5}=[4,6]; L{6}=[4,5];
GSCC = tarjan(L);
assert_mat(length(GSCC),2)
assert_mat(length(GSCC{1}),3)
assert_mat(length(GSCC{2}),3)


for iter=1:100  % completely random matrix testing ....

  
  % undirected graph testing ========================
  adj = [0 1; 0 0];  % initialize so that the while loop does not break
  while not(isConnected(adj)); adj = randomGraph(randi(50)+1,rand); end

  L=adj2adjL(adj);
  GSCC = tarjan(L);
  assert_mat(length(GSCC),1)
  assert_mat(GSCC{1},[1:length(adj)])
  
  % directed graph testing ==========================
  adj=randomDirectedGraph(randi(50)+1,rand);
  L=adj2adjL(adj);
  GSCC = tarjan(L);
  
  
  if isConnected(adj) & isConnected(transpose(adj)) & length(adj)>0
    
    % there should be one component containing all nodes
    assert_mat(length(GSCC),1)
    assert_mat(GSCC{1},[1:length(adj)])
    
    
  else  % disconnected directed graph
    
    ll=[];
    for gg=1:length(GSCC); ll=[ll length(GSCC{gg})]; end;
    [ml,maxll]=max(ll);
    
    assert_mat(isConnected(adj(GSCC{maxll},GSCC{maxll})) | length(GSCC{maxll})==1)
    
    for ii=1:length(adj)
      if isempty(find(GSCC{maxll}==ii))
        
        tryGC = [GSCC{maxll}, ii];
        assert_mat(not(isConnected(adj(tryGC,tryGC))) | not(isConnected(transpose(adj(tryGC,tryGC)))))
        
      end
      
    end
    
  end
  
end
% ================================================
% ================================================


% testing graphComplement.m ======================
fprintf('testing graphComplement.m\n')

mat = [1 0 0 1 1 1; 0 1 0 1 1 1; 0 0 1 0 1 1; 1 1 0 1 0 0; 1 1 1 0 1 0; 1 1 1 0 0 1];
assert_mat(graphComplement(bowtie),mat)
assert_mat(graphComplement(undirected_triangle),eye(3))  
% ================================================


% Testing graphDual.m ============================
fprintf('testing graphDual.m\n')

gd=graphDual(adj2adjL(bowtie));
gdT={};
gdT{1}=[2,3]; gdT{2}=[1,3,4]; gdT{3}=[1,2,4]; gdT{4}=[2,3,5,6]; gdT{5}=[4,6,7]; gdT{6}=[4,5,7]; gdT{7}=[5,6];
assert_mat(gd,gdT)

gd=graphDual(adj2adjL(undirected_triangle));
gdT={};
gdT{1}=[2,3]; gdT{2}=[1,3]; gdT{3}=[1,2];
assert_mat(gd,gdT)

L={}; LT={}; L{1}=[2]; L{2}=[1]; LT{1}=[];
assert_mat(LT,graphDual(L))
% ================================================


% testing subgraph.m =============================
fprintf('testing subgraph.m\n')
assert_mat(undirected_triangle,subgraph(bowtie,[1,2,3]))
% ================================================
  
% testing leafNodes.m ============================
fprintf('testing leafNodes.m\n')
assert_mat(leafNodes(edgeL2adj(undirected_cherry)),[2,3])
assert_mat(leafNodes(edgeL2adj(directed_cherry)),[2,3])
assert_mat(length(leafNodes(undirected_triangle)),0)
% ================================================

% testing leafEdges.m ============================
fprintf('testing leafEdges.m\n')
assert_mat(leafEdges(edgeL2adj(undirected_cherry)),[1,2;1,3])
assert_mat(leafEdges(edgeL2adj(directed_cherry)),[1,2;1,3])
assert_mat(length(leafEdges(undirected_triangle)),0)
hut = [2,1,1;3,1,1];
assert_mat(length(leafEdges(edgeL2adj(hut))),0)
% ================================================


% testing isSimple.m =============================
fprintf('testing isSimple.m\n')

assert_mat(isSimple(randomGraph(randi(5)+20,rand)),true)  % simple graph
assert_mat(isSimple(edgeL2adj([1,2,2])),false)      % multi-edge
assert_mat(isSimple( [1 0 0; 0 0 1; 0 1 0]),false)  % matrix with loops
assert_mat(isSimple([0 1 1; 1 0 0; 0 1 0]),false)   % directed matrix
% ================================================

  
% testing isDirected.m ===========================
fprintf('testing isDirected.m\n')
assert_mat(isDirected(randomDirectedGraph(randi(5)+20,rand)),true)  
assert_mat(isDirected(randomGraph(randi(5)+20,rand)),false)
% ================================================

% testing isSymmetric.m ==========================
fprintf('testing isSymmetric.m\n')

for i=1:100
  assert_mat(isSymmetric(randomGraph(randi(5)+20,rand)),true)

  adj = randomDirectedGraph(randi(5)+20,rand);
  assert_mat(not(isSymmetric(adj)) | adj==zeros(size(adj)) | adj==ones(size(adj)))
end
% ================================================

% testing isConnected.m ==========================
fprintf('testing isConnected.m\n')
assert_mat(isConnected(bowtie),true)
assert_mat(isConnected(disconnected_bowtie),false)
% ================================================

% testing isWeighted.m ===========================
fprintf('testing isWeighted.m\n')
assert_mat(isWeighted([1,2,2]),true)

assert_mat(isWeighted(adj2edgeL(randomGraph(randi(5)+20,rand))),false)
  
assert_mat(isWeighted(adj2edgeL(randomDirectedGraph(randi(5)+20,rand))),false)
  
assert_mat(isWeighted([1,2,0.5; 1,3,1.5; 1,4,1]),true)
assert_mat(isWeighted([1,2,0.5; 1,3,1; 1,4,1]),true)
% ================================================


% testing isRegular.m ============================
fprintf('testing isRegular.m\n')
adj = edgeL2adj(canonicalNets(20,'circle'));
assert_mat(isRegular(adj),true)

adj = edgeL2adj(canonicalNets(20,'tree',3));
assert_mat(isRegular(adj),false)

assert_mat(isRegular([0 1; 1 0]),true)
assert_mat(isRegular([0 0; 1 0]),false)
% ================================================


% testing isComplete.m ===========================
fprintf('testing isComplete.m\n')
assert_mat(isComplete([0 1; 1 0]),true)

assert_mat(isComplete(edgeL2adj(directed_cherry)),false)

assert_mat(isComplete(edgeL2adj(undirected_cherry)),false)

randint = randi(10)+10;
adj = ones(randint)-eye(randint);
assert_mat(isComplete(adj),true)
% ================================================


% testing isEulerian.m ===========================
fprintf('testing isEulerian.m\n')

adj = edgeL2adj(canonicalNets(10,'circle'));
assert_mat(isEulerian(adj),true)

adj = edgeL2adj(canonicalNets(10,'tree',3));
assert_mat(isEulerian(adj),false)
% ================================================

% testing isTree.m ===============================
fprintf('testing isTree.m\n')
adj = edgeL2adj(canonicalNets(randi(10)+10,'tree',2));
assert_mat(isTree(adj),true)

adj = edgeL2adj(canonicalNets(randi(10)+10,'circle'));
assert_mat(isTree(adj),false)
% ================================================

% testing isGraphic.m ============================
fprintf('testing isGraphic.m\n')
for i=1:100
  adj = giantComponent(randomGraph(randi(20)+1,0.5));
  [deg,~,~] = degrees(adj);
  assert_mat(isGraphic(deg) | adj==0)
end
% ================================================


% testing isBipartite.m ==========================
fprintf('testing isBipartite.m\n')

assert_mat(isBipartite(adj2adjL(bowtie)),false)
assert_mat(isBipartite(edgeL2adjL(undirected_cherry)),true)

even_circle = canonicalNets(2*randi(10),'circle');
assert_mat(isBipartite(edgeL2adjL(even_circle)),true)

odd_circle = canonicalNets(2*randi(10)+1,'circle');
assert_mat(isBipartite(edgeL2adjL(odd_circle)),false)
% ================================================


% testing adj2adjL.m =============================
fprintf('testing adj2adjL.m\n')
assert_mat(adj2adjL(bowtie),bowtie_adjL')
% ================================================

% testing adjL2adj.m =============================
fprintf('testing adjL2adj.m\n')

assert_mat(adjL2adj(bowtie_adjL),bowtie)

L = {}; L{1}=[2,3]; L{2}=[]; L{3}=[];
assert_mat(adjL2adj(L),edgeL2adj(directed_cherry))
% ================================================


% testing adj2edgeL.m ============================
fprintf('testing adj2edgeL.m\n')

assert_mat(sortrows(adj2edgeL(bowtie)),bowtie_edgeL)
assert_mat(adj2edgeL([0 1 1; 0 0 0; 0 0 0]),directed_cherry)
% ================================================


% testing edgeL2adj.m ============================
fprintf('testing edgeL2adj.m\n')

assert_mat(edgeL2adj(bowtie_edgeL),bowtie)
assert_mat(edgeL2adj(directed_cherry),[0 1 1; 0 0 0; 0 0 0])
% ================================================


% testing adj2inc.m ==============================
fprintf('testing adj2inc.m\n')

randint = randi(10)+1;
assert_mat(adj2inc(eye(randint)),eye(randint))

assert_mat(adj2inc([0 1 0; 0 1 0; 1 0 0 ]),[-1 0 1; 1 1 0; 0 0 -1])
assert_mat(adj2inc([0 2; 0 0]),[-1 -1; 1 1])  % double edge
% ================================================


% testing inc2adj.m ==============================
fprintf('testing inc2adj.m\n')

randint = randi(10)+1;
assert_mat(inc2adj(eye(randint))==eye(randint))

adj = ones(3) - eye(3);
assert_mat(inc2adj(adj),adj)

inc = [-1 1; 1 0; 0 -1];  % two edges (1->2, 3->1)
assert_mat(inc2adj(inc)==[0 1 0; 0 0 0; 1 0 0])
% ================================================

% testing adj2str.m ==============================
fprintf('testing adj2str.m\n')

assert_mat(adj2str(ones(3)-eye(3)),'.2.3,.1.3,.1.2,')
assert_mat(adj2str(eye(3)),'.1,.2,.3,')
assert_mat(adj2str([0 2; 0 0]),'.2,,')
% ================================================

% testing str2adj.m ==============================
fprintf('testing str2adj.m\n')

assert_mat(ones(3)-eye(3),str2adj('.2.3,.1.3,.1.2,'))
assert_mat(eye(3),str2adj('.1,.2,.3,'))
assert_mat([0 1 0; 0 0 0; 1 0 0 ],str2adj('.2,,.1,'))
% ================================================


% testing adjL2edgeL.m ===========================
fprintf('testing adjL2edgeL.m\n')

assert_mat(adjL2edgeL({[2,3],[],[]}),directed_cherry)
assert_mat(sortrows(adjL2edgeL(bowtie_adjL)),bowtie_edgeL)
% ================================================


% testing edgeL2adjL.m ===========================
fprintf('testing edgeL2adjL.m\n')
assert_mat(edgeL2adjL(directed_cherry),{[2,3],[],[]}')
% ================================================


% testing inc2edgeL.m ============================
fprintf('testing inc2edgeL.m\n')

assert_mat(inc2edgeL([1 0 0; 0 1 0; 0 0 1]),[1 1 1; 2 2 1; 3 3 1])  % three self-loops
assert_mat(inc2edgeL([-1 -1; 1 0; 0 1]),[1 2 1; 1 3 1])
assert_mat(inc2edgeL([-1;1]),[1 2 1])
% ================================================


% testing adj2simple.m ===========================
fprintf('testing adj2simple.m\n')

assert_mat(adj2simple(rand(6)),ones(6)-eye(6))
assert_mat(adj2simple([0 2 0; 1 0 0; 1 2 0]),[0 1 1; 1 0 1; 1 1 0])
% ================================================


% testing edgeL2simple.m =========================
fprintf('testing edgeL2simple.m\n')

assert_mat(length(edgeL2simple([1 1 1; 2 2 1; 3 3 1])),0)
assert_mat(sortrows(edgeL2simple([1 2 1; 1 3 2;4 5 1.4])),[1 2 1; 1 3 1; 2 1 1; 3 1 1; 4 5 1; 5 4 1])
% ================================================


% testing symmetrize.m ===========================
fprintf('testing symmetrize.m\n')
for i=1:20
  adj = randomDirectedGraph(randi(10)+3,rand);
  assert_mat(isSymmetric(symmetrize(adj)),true)
end


% testing symmetrizeEdgeL.m ======================
fprintf('testing symmetrizeEdgeL.m\n')

for x=1:50
  adj = randomDirectedGraph(randi(20)+2,rand); % create a random adjacency
  el = adj2edgeL(adj);
  if isempty(el); continue; end
  elsym = symmetrizeEdgeL(el);
  adjsym = edgeL2adj(elsym);
  assert_mat(isSymmetric(adjsym),true)
end
% ================================================


% testing addEdgeWeights.m =======================
fprintf('testing addEdgeWeights.m\n')

assert_mat([1 2 2; 1 3 1; 3 4 3],addEdgeWeights([1 2 1; 1 2 1; 1 3 1; 3 4 2; 3 4 1]))
assert_mat([1 2 2; 2 3 4],addEdgeWeights([1 2 2; 2 3 4]))
% ================================================


% testing degrees.m ==============================
fprintf('testing degrees.m\n')

assert_mat([2 2 3 3 2 2],degrees(bowtie))
assert_mat([2 1 1],degrees(edgeL2adj(directed_cherry)))
assert_mat([2 1 1],degrees(edgeL2adj(undirected_cherry)))

[deg,indeg,outdeg]=degrees(edgeL2adj(directed_cherry));
assert_mat(deg,[2 1 1])
assert_mat(indeg,[0 1 1])
assert_mat(outdeg,[2 0 0])

assert_mat([4 4 4],degrees([0 2 1; 0 0 1; 1 1 0]))
% ================================================


% testing rewire.m ===============================
fprintf('testing rewire.m\n')

for x=1:100
  
  el = adj2edgeL(randomGraph(randi(10)+10,0.4));
  deg = degrees(edgeL2adj(el));
  eln = rewire(el,randi(5));
  degn = degrees(edgeL2adj(eln));
  
  assert_mat(deg,degn)

end
% ================================================


% testing rewireThisEdge.m =======================
fprintf('testing rewireThisEdge.m\n')

for x=1:100
  
  adj = [0 1; 0 0];
  while not(isConnected(adj)); adj = randomGraph(randi(10)+10,rand); end
  el = adj2edgeL(adj);
  deg = degrees(edgeL2adj(el));
  
  edgeind = randi([1,length(el)]);
  eln = rewireThisEdge(el,el(edgeind,1),el(edgeind,2));
  if isempty(eln); continue; end
  
  adjn = edgeL2adj(eln);
  degn = degrees(adjn);
  
  assert_mat(deg,degn)
  assert_mat(isSimple(adjn),true)

 
end
% ================================================



% testing rewireAssort.m =========================
fprintf('testing rewireAssort.m\n')

for x=1:100
  adj = [0 0; 0 0];
  
  while not(isConnected(adj)); adj = randomGraph(randi(10)+10,0.4); end
  el = adj2edgeL(adj);
  eln = rewireAssort(el,randi(5));
  
  assert_mat(pearson(edgeL2adj(eln))>=pearson(edgeL2adj(el))-10^(-7))
end
% ================================================


% testing rewireDisassort.m ======================
fprintf('testing rewireDisassort.m\n')
for x=1:100
  
  adj = [0 0; 0 0];
  while not(isConnected(adj)); adj = randomGraph(randi(10)+10,0.4); end
  el = adj2edgeL(adj);
  eln = rewireDisassort(el,randi(5));

  assert_mat(pearson(edgeL2adj(eln))<=pearson(edgeL2adj(el))+10^(-7))
  
end
% ================================================

% testing aveNeighborDeg.m =======================
fprintf('testing aveNeighborDeg.m\n')
assert_mat(aveNeighborDeg(undirected_triangle),[2 2 2])
assert_mat(aveNeighborDeg(bowtie),[2.5 2.5 7/3 7/3 2.5 2.5])
% ================================================

% testing sortNodesBySumNeighborDegrees.m ===
fprintf('testing sortNodesBySumNeighborDegrees.m\n')

assert_mat(sortNodesBySumNeighborDegrees(bowtie),[4,3,6,5,2,1]')  
assert_mat(sortNodesBySumNeighborDegrees([0 1 1; 1 0 0; 1 0 0]),[1, 3, 2]')
% ================================================

% testing sortNodesByMaxNeighborDegree.m ====
fprintf('testing sortNodesByMaxNeighborDegree.m\n')

assert_mat(sortNodesByMaxNeighborDegree(bowtie),[4,3,6,5,2,1]')
assert_mat(sortNodesByMaxNeighborDegree(edgeL2adj(undirected_cherry)),[1,3,2]')
% ================================================


% testing closeness.m ============================
fprintf('testing closeness.m\n')
assert_mat(closeness(bowtie)',[1/(1+1+2+3+3), 1/(1+1+2+3+3), 1/(1+1+1+2+2), 1/(1+1+1+2+2), 1/(1+1+2+3+3), 1/(1+1+2+3+3)])
assert_mat(closeness([0 1 1; 1 0 0; 1 0 0]),[0.5 1/3 1/3]')
% ================================================


% testing nodeBetweennessSlow.m ==================
fprintf('testing nodeBetweennessSlow.m\n')
assert_mat(nodeBetweennessSlow([0 1; 1 0]),[0 0])
assert_mat(nodeBetweennessSlow([1 1; 0 0]),[0 0])
assert_mat(nodeBetweennessSlow([0 1 1; 1 0 0; 1 0 0]),[1/3 0 0])
assert_mat(nodeBetweennessSlow(bowtie),[0 0 0.4 0.4 0 0])
% need to test an even cycle eventually: when that works!
% bw = nodeBetweennessSlow(edgeL2adj(canonicalNets(2*randi(10)+2,'circle')));
% assert_mat(bw(1)*ones(1,length(bw)),bw)
% ================================================


% testing nodeBetweennessFaster.m ================
fprintf('testing nodeBetweennessFaster.m\n')
assert_mat(nodeBetweennessFaster([0 1; 1 0]),[0 0])
assert_mat(nodeBetweennessFaster([0 1 1; 1 0 0; 1 0 0]),[1/3 0 0])
assert_mat(nodeBetweennessFaster(bowtie),[0 0 0.4 0.4 0 0])

adj = [0 0; 0 0];
for i=1:100
  
  while not(isConnected(adj)); adj = randomGraph(randi(10)+5,rand); end
  assert_mat(nodeBetweennessSlow(adj),nodeBetweennessFaster(adj))
  
end
% need to test an even cycle eventually: when that works!
% bw = nodeBetweennessFaster(edgeL2adj(canonicalNets(2*randi(10)+2,'circle')));
% assert_mat(bw(1)*ones(1,length(bw)),bw)
% ================================================

% testing edgeBetweenness.m ======================
fprintf('testing edgeBetweenness.m\n')

eb_bowtie = adj2edgeL(bowtie);
eb_bowtie(:,3) = [1/30; 4/30; 1/30; 4/30; 4/30; 4/30; 9/30; 9/30; 4/30; 4/30; 4/30; 1/30; 4/30; 1/30];

assert_mat(edgeBetweenness(bowtie),eb_bowtie)
assert_mat(edgeBetweenness(undirected_triangle),[2 1 1/6; 3 1 1/6; 1 2 1/6; 3 2 1/6; 1 3 1/6; 2 3 1/6])
assert_mat(edgeBetweenness([0 1 1 0; 1 0 1 0; 1 1 0 1; 0 0 1 0]),[2 1 1/12; 3 1 1/6; 1 2 1/12; 3 2 1/6; 1 3 1/6; 2 3 1/6; 4 3 3/12; 3 4 3/12])
% ================================================


% testing eigenCentrality.m ======================
fprintf('testing eigenCentrality.m\n')
[v,~]=eig([0 1 1; 1 0 1; 1 1 0]);
assert_mat(eigenCentrality([0 1 1; 1 0 1; 1 1 0]),v(:,3))
% ================================================


% testing clustCoeff.m ===========================
fprintf('testing clustCoeff.m\n')
assert_mat(clustCoeff(undirected_triangle),1)
assert_mat(clustCoeff(edgeL2adj(undirected_cherry)),0)
assert_mat(clustCoeff(edgeL2adj(canonicalNets(randi(10)+5,'tree',2))),0)
[C1,C2] = clustCoeff(bowtie);
assert_mat([C1,C2],[3/5,7/9])
% ================================================


% testing weightedClustCoeff.m ===================
fprintf('testing weightedClustCoeff.m\n')
randint = randi(20);
assert_mat(length(weightedClustCoeff(randomGraph(randint+5,rand))),randint+5)
% ================================================


% testing pearson.m ==============================
fprintf('testing pearson.m\n')
assert_mat(pearson(star),-1)
assert_mat( pearson( edgeL2adj( canonicalNets(randi(5)+5,'star') ) ) ,-1 )
% ================================================


% testing richClubMetric.m =======================
fprintf('testing richClubMetric.m\n')
assert_mat(richClubMetric(randomGraph(randi(5)+5,rand),12),0)
assert_mat(richClubMetric(bowtie,2),linkDensity(bowtie))

mat = [0 1 1 0; 1 0 1 0; 1 1 0 1; 0 0 1 0];
assert_mat(richClubMetric(mat,2),1)
% ================================================


% testing sMetric.m ==============================
fprintf('testing sMetric.m\n')
assert_mat(sMetric(undirected_triangle),2*12)
assert_mat(sMetric(bowtie),2*41)
assert_mat(sMetric(edgeL2adj(directed_cherry)),4)
assert_mat(sMetric(one_directed_edge),1)
% ================================================


% testing simpleDijkstra.m =======================
fprintf('testing simpleDijkstra.m\n')
assert_mat(simpleDijkstra(bowtie,1),[0, 1, 1, 2, 3, 3])
assert_mat(simpleDijkstra(bowtie,3),[1, 1, 0, 1, 2, 2])

mat = [0 3.5 0 1; 3.5 0 1 0; 0 1 0 1.4; 1 0 1.4 0];
assert_mat(simpleDijkstra(mat,1),[0, 3.4, 2.4, 1])
assert_mat(simpleDijkstra(edgeL2adj(directed_cherry),1),[0, 1, 1])
assert_mat(simpleDijkstra(edgeL2adj(directed_cherry),2),[inf, 0, inf])
% ================================================

% testing dijkstra.m =============================
fprintf('testing dijkstra.m\n')
[d,p]=dijkstra(bowtie,1,5);
assert_mat(d,3)
assert_mat(p,[1,3,4,5])

[d,p]=dijkstra(undirected_triangle,3,[]);
assert_mat(d,[1,1,0])
assert_mat(p,{[3,1],[3,2],[3]})

[d,p] = dijkstra(square,3,[]);
assert_mat(d,[2,1,0,1]);
assert_mat(p,{[3,2,1],[3,2],[3],[3,4]})
% ================================================


% testing shortestPathDP.m =======================
fprintf('testing shortestPathDP.m\n')

[Jb,rb,J,r]=shortestPathDP(bowtie,1,3,size(bowtie,1));
assert_mat(Jb,1)
assert_mat(rb,[1,3])

[Jb,rb,J,r]=shortestPathDP(bowtie,1,4,size(bowtie,1));
assert_mat(Jb,2)
assert_mat(rb,[1,3,4])

[Jb,rb,J,r]=shortestPathDP(bowtie,1,5,size(bowtie,1));
assert_mat(Jb,3)
assert_mat(rb,[1,3,4,5])

[Jb,rb,J,r]=shortestPathDP(edgeL2adj(directed_cherry),1,2,3);
assert_mat(Jb,1)
assert_mat(rb,[1,2])

[Jb,rb,J,r]=shortestPathDP(edgeL2adj(directed_cherry),2,3,3);
assert_mat(Jb,inf)
% ================================================


% test minSpanTree.m =============================
fprintf('testing minSpanTree.m\n')
for x=1:100
  
  adj = [0 1; 0 0];
  while not(isConnected(adj)); adj = randomGraph(randi(50)+5,rand); end
  
  tr = minSpanTree(adj);
  assert_mat(isTree(tr),true)
  assert_mat(length(tr),length(adj));  % tree should have the same
                                   % number of nodes as adj
  
end
% ================================================

% test BFS.m =====================================
fprintf('testing BFS.m\n')

for x=1:100
  adj = [0 1; 0 0];
  while not(isConnected(adj)); adj = randomGraph(randi(50)+5,rand); end

  tr = BFS(adj2adjL(adj),randi(length(adj)));
  tr = symmetrize(adjL2adj(tr));
  assert_mat(isTree(tr),true)
  assert_mat(length(tr),length(adj))
  
end
% ================================================


% testing kneighbors.m ===========================
fprintf('testing kneighbors.m\n')

assert_mat(kneighbors(bowtie,1,3),[1 2 3 4 5 6])
assert_mat(kneighbors(bowtie,3,1),[1 2 4])
assert_mat(kneighbors(undirected_triangle,1,2),[1,2,3])
% ================================================


% testing kminNeighbors.m ========================
fprintf('testing kminNeighbors.m\n')

assert_mat(kminNeighbors(bowtie,1,3),[5, 6])
assert_mat(kminNeighbors(bowtie,3,1),[1, 2, 4])
assert_mat(kminNeighbors(bowtie,3,2),[5, 6])
% ================================================


% testing diameter.m =============================
fprintf('testing diameter.m\n')

assert_mat(diameter(undirected_triangle),1)
assert_mat(diameter(bowtie),3)

el=canonicalNets(randi(10)+5,'line');
adj = edgeL2adj(el);
assert_mat(diameter(adj),length(adj)-1)

el=canonicalNets(randi(10)+5,'circle');
adj = edgeL2adj(el);
assert_mat(diameter(adj),floor(length(adj)/2))
% ================================================


% testing avePathLength.m ========================
fprintf('testing avePathLength.m\n')

assert_mat(avePathLength(bowtie),(0+1+1+2+3+3 +0+1+2+3+3+ 0+1+2+2 +0+1+1 +0+1 +0)/15)
assert_mat(avePathLength(undirected_triangle),1)
adj = edgeL2adj(canonicalNets(6,'line'));
assert_mat(avePathLength(adj),(0+1+2+3+4+5 +0+1+2+3+4 +0+1+2+3 +0+1+2+ 0+1 +0)/15)
% ================================================


% testing smoothDiameter.m =======================
fprintf('testing smoothDiameter.m\n')

adj = [0 1; 0 0];
while not(isConnected(adj)); adj = randomGraph(randi(10)+10,rand); end
assert_mat(diameter(adj),smoothDiameter(adj,1))  % should be the same when the fraction is 1
% ================================================


% testing vertexEccentricity.m ===================
fprintf('testing vertexEccentricity.m\n')

assert_mat(vertexEccentricity(bowtie),[3,3,2,2,3,3])
assert_mat(vertexEccentricity(undirected_triangle),[1,1,1])
% ================================================


% testing graphRadius.m ==========================
fprintf('testing graphRadius.m\n')

assert_mat(graphRadius(bowtie),2)

el = canonicalNets(randi(10)+10,'line');
adj = edgeL2adj(el);
assert_mat(graphRadius(adj),(size(adj,1)-mod(size(adj,1),2))/2)
% ================================================


% testing distanceDistribution.m =================
fprintf('testing distanceDistribution.m\n')

assert_mat(distanceDistribution(bowtie),[7/15, 4/15, 4/15, 0, 0])
assert_mat(distanceDistribution(undirected_triangle),[1, 0])
assert_mat(distanceDistribution(edgeL2adj(undirected_cherry)),[2/3,1/3])
% ================================================


% testing numConnTriples.m =======================
fprintf('testing numConnTriples.m\n')
assert_mat(numConnTriples(bowtie),6)
assert_mat(numConnTriples(undirected_triangle),1)
assert_mat(numConnTriples(edgeL2adj(undirected_cherry)),1)
% ================================================

% testing numLoops.m =============================
fprintf('testing numLoops.m\n')
assert_mat(numLoops(undirected_triangle),1)
assert_mat(numLoops(bowtie),2)
assert_mat(numLoops(edgeL2adj(undirected_cherry)),0)
assert_mat(numLoops(square),1)
% ================================================

% testing loops3.m ===============================
fprintf('testing loops3.m\n')
assert_mat(loops3(bowtie),2)
assert_mat(loops3(square),0)
assert_mat(loops3(undirected_triangle),1)
assert_mat(loops3(edgeL2adj(canonicalNets(randi(10)+3,'btree'))),0)
assert_mat(loops3(edgeL2adj(canonicalNets(4,'trilattice'))),2)
% ================================================


% testing loops4.m ===============================
fprintf('testing loops4.m\n')

assert_mat(loops4(bowtie),{})
c4 = ones(4)-eye(4); % clique of size 4
assert_mat(loops4(c4),{'1-2-3-4'})
c6 = ones(6)-eye(6); % clique of size 6
assert_mat(length(loops4(c6)),nchoosek(6,4))
% ================================================


% testing numStarMotifs.m ========================
fprintf('testing numStarMotifs.m\n')

assert_mat(numStarMotifs(bowtie_adjL,3),4+6)
assert_mat(numStarMotifs(bowtie_adjL,4),2)
assert_mat(numStarMotifs(bowtie_adjL,5),0)

assert_mat(numStarMotifs(adj2adjL(undirected_triangle),3),3)
assert_mat(numStarMotifs(adj2adjL(undirected_triangle),2),6)

assert_mat(numStarMotifs(bowtie_adjL,1),6)   % trivial case
% ================================================


% testing laplacianMatrix.m ======================
fprintf('testing laplacianMatrix.m\n')

assert_mat(laplacianMatrix(bowtie),[2 -1 -1 0 0 0; -1 2 -1 0 0 0; -1 -1 3 -1 0 0; 0 0 -1 3 -1 -1; 0 0 0 -1 2 -1; 0 0 0 -1 -1 2])
assert_mat(laplacianMatrix(undirected_triangle),[2 -1 -1; -1 2 -1; -1 -1 2])
% ================================================


% testing graphSpectrum.m ========================
fprintf('testing graphSpectrum.m\n')
adj = randomGraph(randi(50)+10,rand);
assert_mat(length(graphSpectrum(adj)),length(adj))
% ================================================

% testing algebraicConnectivity.m ================
fprintf('testing algebraicConnectivity.m\n')
adj = randomGraph(randi(50)+10,rand);
assert_mat(length(algebraicConnectivity(adj)),1)
% ================================================

% testing fiedlerVector.m ========================
fprintf('testing fiedlerVector.m\n')
adj = randomGraph(randi(50)+10,rand);
assert_mat(length(fiedlerVector(adj)),length(adj))
% ================================================

% testing graphEnergy.m ==========================
fprintf('testing graphEnergy.m\n')
adj = randomGraph(randi(50)+10,rand);
assert_mat(length(graphEnergy(adj)),1)
% ================================================


% testing simpleSpectralPartitioning.m ===========
fprintf('testing simpleSpectralPartitioning.m\n')

for xx=1:50  % do the randomized test 50 times
  n = randi(99)+11;   % number of nodes
  adj = randomModularGraph(n,4,0.1,3);  % random graph with n nodes
  num_groups = randi(10)+1;  % number of groups to split the nodes in
  groups = [];
  for x=1:length(num_groups)-1; groups = [groups ceil(rand*n/num_groups)+1]; end
  groups = [groups n-sum(groups)];

  modules = simpleSpectralPartitioning(adj,groups);
  for m=1:length(modules); assert_mat(length(modules{m}),groups(m)); end
  
end % end of 50 iterations
% ================================================


% testing newmanGirvan.m =========================
fprintf('testing newmanGirvan.m\n')
modules = newmanGirvan(bowtie,2);
assert_mat(modules{1}==[1,2,3])
assert_mat(modules{2}==[4,5,6])
% ================================================


% testing newmanEigenvectorMethod.m ==============
fprintf('testing newmanEigenvectorMethod.m\n')

modules = newmanEigenvectorMethod(bowtie);
assert_mat(length(modules),2)
assert_mat(modules{1},[4,5,6])
assert_mat(modules{2},[1,2,3])


for x=1:100
  adj = randomGraph(randi(10)+5,1);
  Adj = zeros(4*length(adj));
  Adj(1:length(adj),1:length(adj))=adj;
  Adj(length(adj)+1:2*length(adj),length(adj)+1:2*length(adj))=adj;
  Adj(2*length(adj)+1:3*length(adj),2*length(adj)+1:3*length(adj))=adj;
  Adj(3*length(adj)+1:4*length(adj),3*length(adj)+1:4*length(adj))=adj;

  Adj(5,length(adj)+5)=1; Adj(length(adj)+5,5)=1; 
  Adj(length(adj)+6,2*length(adj)+6)=1; Adj(2*length(adj)+6,length(adj)+6)=1; 
  Adj(2*length(adj)+7,3*length(adj)+7)=1; Adj(3*length(adj)+7,2*length(adj)+7)=1; 
  Adj(3*length(adj)+1,1)=1; Adj(1,3*length(adj)+1)=1; 

  modules = newmanEigenvectorMethod(Adj);
  assert_mat(length(modules),4)


  prescribed = randi(6)+2;
  
  n = randi(50)+50;
  adj = [];
  while not(isConnected(adj)); adj = randomModularGraph(n,prescribed,0.9*log(n)/n,2+0.3*rand); end
  modules = newmanEigenvectorMethod(adj);
  
  sumnodes = 0;
  for m=1:length(modules); sumnodes = sumnodes + length(modules{m}); end
  assert_mat(sumnodes,n)
  
  for m1=1:length(modules)
    for m2=m1+1:length(modules)
      
      assert_mat(length(intersect(modules{m1},modules{m2})),0)
      
    end
  end
    
end
% ================================================


% testing newmanCommFast.m =======================
fprintf('testing newmanCommFast.m\n')

[gH,Q]=newmanCommFast(bowtie);
close all;
assert_mat(max(Q),Q(6-1));

[gH,Q]=newmanCommFast(randomModularGraph(100,4,0.1,10));
close all;
assert_mat(length(gH),length(Q))
[~,ind]=max(Q);
assert_mat(length(gH{ind}),4)
% ================================================


% testing modularityMetric.m =====================
fprintf('testing modularityMetric.m\n')

for i=1:20
  
  adj = [0 1; 0 0];
  num_modules = randi([2,5]);
  while not(isConnected(adj)); adj = randomModularGraph(30,num_modules,0.1,5); end
 
  % compare to newmanCommFast
  [mH,Q1] = newmanCommFast(adj);
  close all;
  Q2 = [];
  for m=1:length(mH); Q2 = [Q2 modularityMetric(mH{m},adj)]; end
  
  assert_mat(Q1,Q2)

  % compare to the newman-girvan routine
  [modules0,~,Q0] = newmanGirvan(adj,num_modules);
  assert_mat(Q0,modularityMetric(modules0,adj))
 
end
% ================================================


% testing louvainCommunityFinding.m ==============
fprintf('testing louvainCommunityFinding.m\n');

extended_bowtie0 = [1 2; 2 3; 3 2; 3 4; 4 5; 5 6; 4 6; 6 7; 7 8; 7 9; 8 9];
extended_bowtie = [];
for row=1:size(extended_bowtie0,1)
  extended_bowtie = [extended_bowtie; extended_bowtie0(row,:) 1];
end
clear extended_bowtie0
extended_bowtie = symmetrizeEdgeL(extended_bowtie);
adj = edgeL2adj(extended_bowtie);

[modules,inmodule]=louvainCommunityFinding(adj);
assert_mat(length(modules),3)
assert_mat([inmodule{1},inmodule{2},inmodule{3}],[1,1,1]*inmodule{1})
assert_mat([inmodule{4},inmodule{5},inmodule{6}],[1,1,1]*inmodule{4})
assert_mat([inmodule{7},inmodule{8},inmodule{9}],[1,1,1]*inmodule{7})

[modules,inmodule]=louvainCommunityFinding(bowtie);
assert_mat(length(modules),2)
assert_mat([inmodule{1},inmodule{2},inmodule{3}],[1,1,1]*inmodule{1})
assert_mat([inmodule{4},inmodule{5},inmodule{6}],[1,1,1]*inmodule{4})

% concatenate 4 complete graphs: 
adj = ones(10,10)-eye(10);
Adj = zeros(40,40);

Adj(1:10,1:10)=adj;
Adj(11:20,11:20)=adj;
Adj(21:30,21:30)=adj;
Adj(31:40,31:40)=adj;

Adj(10,11) = 1; Adj(11,10) = 1;
Adj(20,21) = 1; Adj(21,20) = 1;
Adj(30,31) = 1; Adj(31,30) = 1;

[modules,inmodule]=louvainCommunityFinding(Adj);
assert_mat(length(modules),4)

% concatenate 4 dense graphs
adj = [0 1; 0 0];
while not(isConnected(adj)); adj = randomGraph(10,0.9); end
Adj = zeros(40,40);

Adj(1:10,1:10)=adj;
Adj(11:20,11:20)=adj;
Adj(21:30,21:30)=adj;
Adj(31:40,31:40)=adj;

Adj(10,11) = 1; Adj(11,10) = 1;
Adj(20,21) = 1; Adj(21,20) = 1;
Adj(30,31) = 1; Adj(31,30) = 1;

[modules,inmodule]=louvainCommunityFinding(Adj);
assert_mat(length(modules),4)
% ================================================



% testing randomGraph.m ==========================
fprintf('testing randomGraph.m\n');

% testing the size of the graph
randint = randi(20)+3;
assert_mat(size(randomGraph(randint),1),randint)
assert_mat(size(randomGraph(randint),2),randint)

% testing the default probability of attachment
for x=1:50
  randint = randi(50)+50;
  adj = randomGraph(randint);
  assert_mat(linkDensity(adj)>0.35);
  assert_mat(linkDensity(adj)<0.65);
end

% testing a random probability of attachment
for x=1:50
  p = rand;
  randint = randi(50)+50;
  adj = randomGraph(randint,p);
  assert_mat(linkDensity(adj)>p-0.15);
  assert_mat(linkDensity(adj)<p+0.15);
end

% testing for the number of edges, E
for x=1:50
  randint = randi(50)+50;
  E = randi([1,randint-1]);
  adj = randomGraph(randint,[],E);
  assert_mat(numEdges(adj),E);
end
% ================================================


% testing randomDirectedGraph.m ==================
fprintf('testing randomDirectedGraph.m\n');

for i=1:30
  p=rand;
  n = randi(40)+40;
  adj = randomDirectedGraph(n,p);
  assert_mat(linkDensity(adj)>p-0.15)
  assert_mat(linkDensity(adj)<p+0.15)

  assert_mat(size(adj),[n,n]);
end
% ================================================


% test graphFromDegreeSequence.m =================
fprintf('testing graphFromDegreeSequence.m\n')

for x=1:50
  adj = [0 1; 0 0];
  while not(isConnected(adj)); adj = randomGraph(randi(50)+50,rand); end
  adjr=graphFromDegreeSequence(degrees(adj));
  assert_mat(isSimple(adjr),true)
  assert_mat(degrees(adj),degrees(adjr))
end
% ================================================


% test randomGraphFromDegreeSequence.m ===========
fprintf('testing randomGraphFromDegreeSequence.m\n')

for x=1:100
  
  adj = [0 1; 0 0];
  N = randi(50)+10;
  while not(isConnected(adj)); adj = randomGraph(N,log(N)/N); end
  
  adjr = randomGraphFromDegreeSequence(degrees(adj));
  
  assert_mat(isSimple(adjr),true)
  assert_mat(degrees(adj),degrees(adjr))
end
% ================================================


% test canonicalNets.m ===========================
fprintf('testing canonicalNets.m\n');

for x=1:20
  N = randi(50)+10;
  
  elLine = canonicalNets(N,'line'); % test line
  adj = edgeL2adj(elLine);
  assert_mat(numNodes(adj),N);
  assert_mat(isConnected(adj),true)
  assert_mat(degrees(adj),[1 2*ones(1,N-2) 1])
  assert_mat(isTree(adj),true)
  
  elC = canonicalNets(N,'circle'); % test circle
  adj = edgeL2adj(elC);
  assert_mat(numNodes(adj),N);
  assert_mat(isConnected(adj),true)
  assert_mat(degrees(adj),2*ones(1,N))
  assert_mat(isTree(adj),false)
  
  elS = canonicalNets(N,'star'); % test star
  adj = edgeL2adj(elS);
  assert_mat(numNodes(adj),N);
  assert_mat(isConnected(adj),true)
  assert_mat(degrees(adj),[N-1 ones(1,N-1)])
  assert_mat(isTree(adj),true)
    
  elCl = canonicalNets(N,'clique'); % test clique
  adj = edgeL2adj(elCl);
  assert_mat(numNodes(adj),N)
  assert_mat(isComplete(adj),true)
  assert_mat(degrees(adj),(N-1)*ones(1,N))
  
  
  elBT = canonicalNets(N,'btree'); % test binary tree
  adj = edgeL2adj(elBT);
  assert_mat(numNodes(adj),N)
  assert_mat(isTree(adj),true)
  deg = degrees(adj);
  assert_mat(max(deg),3)
  assert_mat(deg(1),2)
  assert_mat(deg(N),1)
  
  b = randi([3,6]);
  elT = canonicalNets(N,'tree',b); % test general tree
  adj = edgeL2adj(elT);
  assert_mat(numNodes(adj),N)
  assert_mat(isTree(adj),true)
  deg = degrees(adj);
  assert_mat(max(deg)<=b+1,true)
  assert_mat(deg(1),b)
  assert_mat(deg(N),1)
  
  assert_mat(edgeL2adj(canonicalNets(N,'tree',2)),edgeL2adj(canonicalNets(N,'btree')))
  
  b = randi([3,6]);
  elH = canonicalNets(N,'hierarchy',b); % test hierarchy
  adj = edgeL2adj(elH);  
  assert_mat(numNodes(adj),N)
  assert_mat(isTree(adj),false)
  deg = degrees(adj);
  assert_mat(max(deg)<=b+3,true)
  assert_mat(deg(1),b)
  
  
  eltr = canonicalNets(N,'trilattice');  % test triangular lattice
  adj = edgeL2adj(eltr);
  assert_mat(numNodes(adj),N)
  assert_mat(isTree(adj),false)
  assert_mat(isComplete(adj),false)
  assert_mat(max(degrees(adj))<=6,true)
  assert_mat(isConnected(adj),true)

  
  elsq = canonicalNets(N,'sqlattice');  % test triangular lattice
  adj = edgeL2adj(elsq);
  assert_mat(numNodes(adj),N)
  assert_mat(isTree(adj),false)
  assert_mat(isComplete(adj),false)
  assert_mat(max(degrees(adj))<=4,true)
  assert_mat(isConnected(adj),true)

  elhex = canonicalNets(N,'hexlattice');  % test triangular lattice
  adj = edgeL2adj(elhex);
  assert_mat(numNodes(adj),N)
  assert_mat(isTree(adj),false)
  assert_mat(isComplete(adj),false)
  assert_mat(max(degrees(adj))<=3,true)
  assert_mat(isConnected(adj),true)

  
end
% ================================================


% test kregular.m ================================
fprintf('testing kregular.m\n');
for x=1:30
  
  n = randi(20)+5;   % random integer between 6 and 25
  k = randi(n-2)+1;  % randon integer between 2 and n-1
  if mod(k,2)==1 & mod(n,2)==1; continue; end  % no solution for this case
  el = kregular(n,k);
  adj = edgeL2adj(el);
  assert_mat(degrees(adj),k*ones(1,length(adj)))

end
% ================================================


% test randomGraphDegreeDist.m ===================
fprintf('testing randomGraphDegreeDist.m\n')

N = randi(20)+10;
adj = randomGraphDegreeDist(N,'uniform');
assert_mat(numNodes(adj),N)
assert_mat(isSimple(adj),true)

%adj = randomGraphDegreeDist(N,'normal');
%assert_mat(numNodes(adj),N)
%assert_mat(isSimple(adj),true)

%adj = randomGraphDegreeDist(N,'binomial');
%assert_mat(numNodes(adj),N)
%assert_mat(isSimple(adj),true)

adj = randomGraphDegreeDist(N,'exponential');
assert_mat(numNodes(adj),N)
assert_mat(isSimple(adj),true)

%adj = randomGraphDegreeDist(6,'custom',[1/5 1/5 1/5 1/5 1/5]);
%assert_mat(numNodes(adj),6)
%assert_mat(isSimple(adj),true)

adj = randomGraphDegreeDist(N,'custom');
assert_mat(isempty(adj),true)

adj = randomGraphDegreeDist(N,'anything here');
assert_mat(isempty(adj),true)
% ================================================


% test randomModularGraph.m ======================
fprintf('testing randomModularGraph.m\n');
for x=1:20
  N = randi(50)+10;
  c = randi(5)+1;
  [adj, modules] = randomModularGraph(N,c,0.2,4);
  assert_mat(numNodes(adj),N)
  assert_mat(length(modules),c)
end
% ================================================


% testing weightedRandomSample.m =================
fprintf('testing weightedRandomSample.m\n');

for x=1:20
  n=randi(10)+1;  % number of draws
  P = [1:randi(7)+10];  % population to draw from
  W = rand(1,length(P));
  W = W/sum(W);         % normalize to create weights that sum to 1
  s = weightedRandomSample(n,P,W);
  assert_mat(length(s),n)
  assert_mat(intersect(P,s),unique(s))   % test that all draws exist in the population                                     
end
% ================================================


% testing buildSmaxGraph.m =======================
fprintf('testing buildSmaxGraph.m\n');

for x=1:50
  adj  = [];
  while not(isConnected(adj)); adj = randomGraph(20,0.1); end
  sm = sMetric(adj);
  elmax1 = buildSmaxGraph(degrees(adj));
  adjmax1 = symmetrize(edgeL2adj(elmax1));
  
  smax = sMetric(adjmax1);
  assert_mat(degrees(adjmax1),degrees(adj))
  assert_mat(smax>=sm,true)
  
  elmax2 = buildSmaxGraph(degrees(adj));
  assert_mat(elmax2,elmax1)
end
% ================================================


% testing PriceModel.m ===========================
fprintf('testing PriceModel.m\n')
for x=1:20
  randint = randi(10)+10;
  adj = PriceModel(randint);
  assert_mat(isDirected(adj),true)
  assert_mat(numNodes(adj),randint)    
end
% ================================================


% testing preferentialAttachment.m ===============
fprintf('testing preferentialAttachment.m\n')
for x=1:10
  el = preferentialAttachment(randi(10)+10,1);
  adj = edgeL2adj(el);
  assert_mat(isTree(adj),true)
  assert_mat(isSimple(adj),true)
  
  randint = randi(30)+5;
  el = preferentialAttachment(randint,2);
  adj = edgeL2adj(el);
  assert_mat(numEdges(adj),1+2*(length(adj)-2))
    
end
% ================================================


% testing exponentialGrowthModel.m ===============
fprintf('testing exponentialGrowthModel.m\n')
for x=1:10
  el = exponentialGrowthModel(randi(100));
  adj=edgeL2adj(el);
  assert_mat(isConnected(adj),true)
  assert_mat(isTree(adj),true)
end
% ================================================

% testing masterEquation.m =======================
fprintf('testing masterEquation.m\n')
for x=1:30
  randint = randi(100)+3;
  adj = masterEquationGrowthModel(randint,1,0);
  assert_mat(isTree(adj),true)
  
  adj = masterEquationGrowthModel(randint,2,0);
  assert_mat(isTree(adj),false)
  
  adj = masterEquationGrowthModel(randint,2,2);
  assert_mat(isSimple(adj),true)
  
end
% ================================================


% testing newmanGastner.m ========================
fprintf('testing newmanGastner.m\n')

for x=1:10
  N = randi(100)+10;
  el = newmanGastner(N,rand,[],'off');  % no plot
  adj = symmetrize(edgeL2adj(el));
  assert_mat(numNodes(adj),N);
  assert_mat(isSimple(adj),true)
end
% ================================================


% testing fabrikantModel.m =======================
fprintf('testing fabrikantModel.m\n')
for x=1:20
  adj = fabrikantModel(randi(30)+10,rand*10,'off');
  assert_mat(isConnected(adj),true)
  assert_mat(isTree(adj),true)
  assert_mat(isSimple(adj),true)
end
% ================================================

% testing DoddsWattsSabel.m ======================
fprintf('testing DoddsWattsSabel.m\n')
for x=1:40
  randint = randi(50)+2;
  m = randi(round(randint/4));
  adj = DoddsWattsSabel(randint,2,m,10*rand,10*rand);
  assert_mat(numEdges(adj),m+randint-1)
  assert_mat(isTree(adj),false)
  assert_mat(isConnected(adj),true)
  assert_mat(isSimple(adj),true)
end
% ================================================

% testing nestedHierarchiesModel.m ===============
fprintf('testing nestedHierarchiesModel.m\n')

el = nestedHierarchiesModel(640,3,[10, 20, 40],10);
adj = edgeL2adj(el);
assert_mat(isSimple(adj));
% ================================================


% testing forestFireModel.m ======================
fprintf('testing forestFireModel.m\n');

for x=1:20
  randint = randi(20)+5;
  L = forestFireModel(randint,rand,10*rand);
  adj = symmetrize(adjL2adj(L));
  assert_mat(isSimple(adj),true)
  assert_mat(randint,numNodes(L))

end
% ================================================

% testing pdfCdfRank.m ===========================
fprintf('testing pdfCdfRank.m\n');
adj = randomGraph(randi(30)+30,0.2);
[xp,yp,xc,yc,lk,lx] = pdfCdfRank(degrees(adj),'off');
assert_mat(length(xp),length(xc))
assert_mat(length(xp),length(yp))
assert_mat(length(yp),length(yc))
assert_mat(length(lk),length(lx))
% ================================================
