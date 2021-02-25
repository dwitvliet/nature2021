function C = consensus_consistency(agreeMat)
% consistency of agreement matrix
%https://an.kaist.ac.kr/~haewoon/papers/2009-imc-consistent.pdf
% eq. (4)

numE = sum(sum(agreeMat>0)) ;
E = agreeMat(agreeMat>0) ; 
C = sum((E - 0.5).^2) / numE * 1/0.25 ;
