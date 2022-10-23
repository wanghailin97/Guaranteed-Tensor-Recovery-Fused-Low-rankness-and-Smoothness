function Zvals = cp_wfg_sparse_setup(Z,W)
%CP_WFG_SPARSE_SETUP Creates a special array.
%
%   ZVALS = CP_WFG_SPARSE_SETUP(Z,W) creates an array ZVALS that
%   contains the values of Z corresponding to the indices specified
%   by W.subs. 
%
%   See also CP_WFG_SPARSE.
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
% $Id: cp_wfg_sparse_setup.m,v 1.2 2010/03/19 23:46:32 tgkolda Exp $

Zsubs = Z.subs;
Wsubs = W.subs;
Zvals = zeros(size(W.vals));
[junk,loc] = ismember(Zsubs,Wsubs,'rows');
Zvals(loc) = Z.vals;
