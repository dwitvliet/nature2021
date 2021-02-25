function [ randizeWSBMModel ] = randomize_wsbm_para(origModel,randParam)
% randomize the parameters of WSBM
% randParam 1=edge, 2=weight, 3=both, random permute

if nargin < 2
    disp('two args')
    return
end

e_perm_idx = randperm(size(origModel.Para.theta_e,1)) ;
w_perm_idx = randperm(size(origModel.Para.theta_w,1)) ;

switch randParam
    case 1     
        new_theta_e = origModel.Para.theta_e(e_perm_idx,:) ;
        new_theta_w = origModel.Para.theta_w ;
    case 2
        new_theta_e = origModel.Para.theta_e ;
        new_theta_w = origModel.Para.theta_w(w_perm_idx,:) ;
    case 3        
        new_theta_e = origModel.Para.theta_e(e_perm_idx,:) ;
        new_theta_w = origModel.Para.theta_w(w_perm_idx,:) ;    
end

randizeWSBMModel = struct(origModel) ; 
randizeWSBMModel.Para.theta_e = new_theta_e;
randizeWSBMModel.Para.theta_w = new_theta_w;







