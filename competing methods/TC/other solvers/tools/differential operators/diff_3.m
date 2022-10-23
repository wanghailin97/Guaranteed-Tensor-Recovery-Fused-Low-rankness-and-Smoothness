function [ DX_3 ] = diff_3(X)
%DIFF_3 
slice = X(:,:,1)- X(:,:,end);
DX_3 = diff(X,1,3);
DX_3 = cat(3,DX_3,slice);

