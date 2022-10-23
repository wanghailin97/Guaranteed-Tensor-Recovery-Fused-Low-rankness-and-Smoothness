function [X] = Unfold( X, dim, i )
%dim是一个维度的向量，i是要沿第几个维展开
X = reshape(shiftdim(X,i-1), dim(i), []);