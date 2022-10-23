function [ DX_1T ] = diff_1T(X)
%DIFF_1T
slice = X(1,:,:)- X(end,:,:);
DX_1T = diff(X,1,1);
DX_1T = cat(1,-slice,-DX_1T);



