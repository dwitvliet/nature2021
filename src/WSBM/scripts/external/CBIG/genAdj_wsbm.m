function [genMat,noNaGenMat] = genAdj_wsbm(wsbmModel)
% generate a synthetic mat given the parameters of a already realized wsbm
% model fit
%
% function [Edge_List,True_Model] = generateEdges(W_truth,E_truth,R,theta_w,theta_e,group_Size,degree_Para)

weightDist_model = wsbmModel.W_Distr ;
edgeDist_model = wsbmModel.E_Distr ;
r_struct_model = wsbmModel.R_Struct.R ;

theta_w_model = wsbmModel.Para.theta_w ;
% tmp = triu(theta_w_model,1) ;
% theta_w_model = theta_w_model + tmp ;

theta_e_model = wsbmModel.Para.theta_e ;
% tmp = triu(theta_w_model,1) ;
% theta_w_model = theta_w_model + tmp ;

[~,tmp] = make_WSBM_prior(wsbmModel) ;
groupS_model = sum(tmp,2) ;

% edgeList = generateEdges_nonNeg(...
%     weightDist_model, edgeDist_model, r_struct_model, ...
%     theta_w_model, theta_e_model, groupS_model) ;

edgeList = generateEdges(...
    weightDist_model, edgeDist_model, r_struct_model, ...
    theta_w_model, theta_e_model, groupS_model) ;

genMat = Edg2Adj(edgeList) ;

% symmetrize the genMat
genMat = triu(genMat) + triu(genMat,1)';

if nargout > 1
    noNaGenMat =  genMat ;
    noNaGenMat(~~isnan(genMat)) = 0 ;
end








