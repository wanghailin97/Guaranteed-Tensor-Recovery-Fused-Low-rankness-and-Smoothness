function s = spmatrix(a)
%SPMATRIX Converts a two-way sparse tensor to sparse matrix.
%
%   SPMATRIX(X) converts a sparse tensor to a sparse matrix. The sparse
%   tensor must be two-dimensional.
%
%   See also SPTENSOR, SPTENSOR/RESHAPE, SPTENMAT
%
%MATLAB Tensor Toolbox.
%Copyright 2010, Sandia Corporation. 

% This is the MATLAB Tensor Toolbox by Brett Bader and Tamara Kolda. 
% http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2010) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in tensor_toolbox/LICENSE.txt
% $Id: spmatrix.m,v 1.4 2010/03/19 23:46:30 tgkolda Exp $

if ndims(a) ~= 2
    error('Sparse tensor must be two dimensional.');
end


if isempty(a.subs)
    s = sparse(a.size(1), a.size(2));
else
    s = sparse(a.subs(:,1), a.subs(:,2), a.vals, a.size(1), a.size(2));
end
