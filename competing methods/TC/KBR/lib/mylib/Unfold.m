function [X] = Unfold( X, dim, i )
%dim��һ��ά�ȵ�������i��Ҫ�صڼ���άչ��
X = reshape(shiftdim(X,i-1), dim(i), []);