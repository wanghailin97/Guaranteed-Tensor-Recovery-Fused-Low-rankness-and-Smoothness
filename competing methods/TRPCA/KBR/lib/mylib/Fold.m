function [X] = Fold(X, dim, i)
%dim��һ��ά�ȵ�������i��Ҫ�صڼ���άչ��
dim = circshift(dim, [1-i, 1-i]);
X = shiftdim(reshape(X, dim), length(dim)+1-i);