function [f,G] = cp_wfg(Z,W,A,normZsqr)
%CP_WFG Function and gradient of CP with missing data.
%
%   [F,G] = CP_WFG(Z,W,A) computes the function and gradient values of
%   the function 0.5 * || W .* (Z - ktensor(A)) ||^2. The input A is a
%   cell array containing the factor matrices. The input W is a (dense
%   or sparse) tensor containing zeros wherever data is missing. The
%   input Z is a (dense or sparse) tensor that is assumed to have
%   zeros wherever there is missing data. The output is the function F
%   and a cell array G containing the partial derivatives with respect
%   to the factor matrices.
%
%   [F,G] = CP_WFG(Z,W,A,NORMZSQR) also passes in the pre-computed
%   norm of Z, which makes the computations faster. 
%
%   See also CP_WFUN, CP_WFG_SPARSE, CP_WFG_SPARSE_SETUP.
%
%MATLAB Tensor Toolbox.
%Copyright 2010, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by Brett Bader, Tamara Kolda,
% Evrim Acar, and Daniel Dunlavy.
% http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2010) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in tensor_toolbox/LICENSE.txt
% $Id: cp_wfg.m,v 1.2 2010/03/19 23:46:32 tgkolda Exp $

%% Compute B = W.*ktensor(A)
if isa(W,'sptensor')
    B = W.*ktensor(A);
else
    B = W.*full(ktensor(A));
end

%% Compute normZ
if ~exist('normZsqr','var')
    normZsqr = norm(Z)^2;
end

% function value
f = 0.5 * normZsqr - innerprod(Z,B) + 0.5 * norm(B)^2;

% gradient computation
N = ndims(Z);
G = cell(N,1);
T = Z - B;
for n = 1:N
    G{n} = zeros(size(A{n}));
    G{n} = -mttkrp(T,A,n);
end

