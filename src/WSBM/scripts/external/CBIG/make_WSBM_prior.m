function [ prior , harshPrior ] = make_WSBM_prior(iMu,wieghtMore)
% takes in a WSBM model, and ouputs at prior matrix based of the mu_0
% of the WSBM model provided, and a harsh prior with just the 
% 'winners' of the mu_0 mat

% lets make the prior assignment matrix
% initialize the matrix
%priorMat = ones(wsbmModel.R_Struct.k , length(wsbmModel.Para.mu)) ;

if nargin < 2
    wieghtMore = 0 ; 
end

% if inputMu is a struct, extract relevant variables
if isstruct(iMu)
    k = size(iMu.Para.mu,1);
    mu = iMu.Para.mu ;
else
    k = size(iMu,1) ;
    mu = iMu ; 
end

%twice as likely (1 would be 100%, double more likely)

weight_label = (1 + wieghtMore ) / (k + wieghtMore) ;
weight_other = 1 / ( k + wieghtMore ) ; 

% % and then fill the mu_prior initially with the weights_others
prior = (ones( k , length(mu))) * weight_other ;

%sanity check 
if (weight_label < weight_other) 
    disp('prior weight is too low dude')
    return
end

% [~,ind] = max(mu,[],1) ;
% harshPrior = dummyvar(ind)' ;
% 
% prior = harshPrior .* 1 ;
% prior(prior==1) = weight_label ;
% prior(prior==0) = weight_other ;

%loop through the best_model.para.mu 
for i = 1:length(mu),
    [~,pos] = max(mu(:,i)) ; 
    prior(pos,i) = weight_label ;
end

%lets alos make the harsh prior
harshPrior = zeros( k , length(mu)) ;
for i = 1:length(mu),
    [~,pos] = max(mu(:,i)) ; 
    harshPrior(pos,i) = 1 ;
end


