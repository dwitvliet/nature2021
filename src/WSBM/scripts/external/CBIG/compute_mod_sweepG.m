function [ sweepResults ] = compute_mod_sweepG( inputData , gammaRange )
% compute modularity here...

if nargin < 2
    gammaRange = 0.5:0.01:4.0 ;
end

%% data
% adjacency mat
if isstruct(inputData)
    M = inputData.Data.Raw_Data ;
else
    M = inputData ;
end

M(isnan(M)) = 0;
N = size(M,1);
R = 2;

% prealocate 
ci_r = zeros(N,R);
q_r = zeros(1,R);

lg = length(gammaRange);

%% run stuff

% loop over gamma vals
for g=1:lg
    %disp(num2str(g));
    gam = gammaRange(g);
    for r=1:R
        [ci_r(:,r) , q_r(r)] = modularity_und(M,gam);
    end;
    % get max q
    ff = find(q_r==max(q_r)); 
    ff = ff(1);
    q_g(g) = q_r(ff);
    ci_g(:,g) = ci_r(:,ff);
end;

% % plot number of modules across gamma range
% figure
% plot(max(ci_g));

% % get stability
% VI = zeros(lg,lg); MI = zeros(lg,lg);
% for i=1:lg
%     for j=1:lg
%         [VI(i,j) , MI(i,j)] = partition_distance(ci_g(:,i),ci_g(:,j));
%     end;
% end;

% % plot the VI across all gamma levels - blue patches are "islands of stable
% % solutions" - use this to select interesting levels of gamma
% figure('position',[50 50 800 800]);
% imagesc(VI);

% % compare to SBM output
% [~ , jj] = max(templateModel.Para.mu);
% SBM = jj';  % not sure this is correct - what I want is the block assignment vector
% for i=1:lg
%     [VI_SBM(:,i) , MI_SBM(:,i)] = partition_distance(ci_g(:,i),SBM);
% end;
% % % plot MI of ci(gamma) versus SBM
% % figure
% % plot(VI_SBM);

% %% pick one gamma val to display
% 
% % display matrix
% % pick the gamma val provided... or use default 1, set up top
% g_display = find(gamvals == gammaVal) ;
% if isempty(g_display)
%    disp('bad gamma val provided')
%    return
% end
% CI = ci_g(:,g_display);
% nummod = max(CI);
% [~ , bb] = sort(CI);

% % display modules
% figure('position',[100 50 800 800]);
% imagesc(log(M(bb,bb)));         % log scale!
% mod_boundaries = find(abs(diff(CI(bb)))>0);
% for i=1:nummod-1
%     ll = line([mod_boundaries(i)+0.5 mod_boundaries(i)+0.5],[0.5 N+0.5]); set(ll,'Color',[1 1 1],'LineWidth',1);
%     ll = line([0.5 N+0.5],[mod_boundaries(i)+0.5 mod_boundaries(i)+0.5]); set(ll,'Color',[1 1 1],'LineWidth',1);
% end;
% colormap(jet); axis square;

%% return the modular communities...
% repends on gamma choice for now

% numNodes = size(inputData.Para.mu,2) ;
% communityLabels_mod = [ (1:numNodes)' CI ] ;

% return mod sweep results
% make top row modularity coeffs
sweepResults = vertcat( q_g , ci_g ) ;

