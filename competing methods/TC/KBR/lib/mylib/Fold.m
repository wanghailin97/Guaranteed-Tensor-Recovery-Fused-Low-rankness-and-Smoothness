function [X] = Fold(X, dim, i)
%dim是一个维度的向量，i是要沿第几个维展开
dim = circshift(dim, [1-i, 1-i]);
X = shiftdim(reshape(X, dim), length(dim)+1-i);