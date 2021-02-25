function [ a ] = sym_RStruct(size)
% return the symmetic RStruct based on the size

a = reshape(linspace(1, size^2, size^2), size, size)';