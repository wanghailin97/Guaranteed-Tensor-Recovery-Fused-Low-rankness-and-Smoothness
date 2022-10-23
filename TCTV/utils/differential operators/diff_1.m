function [ DX_1 ] = diff_1(X)
%DIFF_1 
slice = X(1,:,:)- X(end,:,:);
DX_1  = diff(X,1,1);
DX_1  = cat(1,DX_1,slice);

