function [ a ] = diagWBlocks_RStruct(size)
% return the symmetic RStruct based on the size

a = triu(ones(size),1) ;
b = (size * (size - 1) / 2 ) + 1 ;
a(~~a) = 2:b ;
a = a + eye(size) ;
a = a + triu(a,1)' ;