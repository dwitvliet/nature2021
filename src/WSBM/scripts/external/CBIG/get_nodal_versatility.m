function [V] = get_nodal_versatility(comMat)
% com mat should be nodes x paritions

CM = agreement(comMat) ./ size(comMat,2);

Cs = sin(pi*CM);
V = sum(Cs, 1)./sum(CM, 1);
V(V<1e-10) = 0;
    
    