function [ communityAssign, justLabels ] = community_assign(muLabels)
%% make community assign from muLabels

% if wsbm struct provided, get mu labels out of it
if isstruct(muLabels) 
   muLabels = muLabels.Para.mu;
end

communityAssign = zeros(length(muLabels),2) ; 

for i = 1:length(muLabels),
    
    %make first column the label number
    % i dont think the first column even matters...
    communityAssign(i,1) = i ; %select_nodes(i) ; 
    
    %make second column the 
    [~,temp] = max(muLabels(:,i)) ; 
    communityAssign(i,2) = temp ; 
    
end
    
disp('the block assisgnments of Init Model are:\n')
disp(communityAssign)

if nargout == 2
    justLabels = communityAssign(:,2);
end