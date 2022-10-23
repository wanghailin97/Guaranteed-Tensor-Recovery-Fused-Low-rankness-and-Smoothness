function [ DX_2 ] = diff_2(X)
%DIFF_2
slice = X(:,1,:)- X(:,end,:);
DX_2 = diff(X,1,2);
DX_2 = cat(2,DX_2,slice);

