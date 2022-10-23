function [ DX_3T ] = diff_3T(X)
%DIFF_3T 
slice = X(:,:,1)- X(:,:,end);
DX_3T = diff(X,1,3);
DX_3T = cat(3,-slice,-DX_3T);