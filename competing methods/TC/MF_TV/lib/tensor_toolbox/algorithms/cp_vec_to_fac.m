function A = cp_vec_to_fac(x,Z)
%CP_VEC_TO_FAC Converts a vector to a cell array of factor matrices.
%
%   A = CP_VEC_TO_FAC(X,Z) converts the vector X into a cell array
%   of factor matrices consistent with the size of the tensor Z.
%
%   See also FAC_TO_VEC, CP_FUN, CP_OPT.
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
% $Id: cp_vec_to_fac.m,v 1.3 2010/03/19 23:46:32 tgkolda Exp $

%% Set-up
P = length(x);
N = ndims(Z);
sz = size(Z);

%% Determine R
R = P / sum(sz);

%% Create A
A = cell(N,1);
for n = 1:N
    idx1 = sum(sz(1:n-1))*R + 1;
    idx2 = sum(sz(1:n))*R;
    A{n} = reshape(x(idx1:idx2),sz(n),R);
end