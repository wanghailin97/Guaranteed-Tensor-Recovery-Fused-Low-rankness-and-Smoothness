function [ DX_2T ] = diff_2T(X)
%DIFF_2T
slice = X(:,1,:)- X(:,end,:);
DX_2T = diff(X,1,2);
DX_2T = cat(2,-slice,-DX_2T);