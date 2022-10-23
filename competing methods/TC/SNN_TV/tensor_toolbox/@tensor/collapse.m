function Y = collapse(X,dims,fun)
%COLLAPSE Collapse tensor along specified dimensions.
%
%   Y = COLLAPSE(X,DIMS) sums the entries of X along all dimensions
%   specified in DIMS. If DIMS is negative, then X is summed across
%   all dimensions *not* specified by -DIMS.
%
%   Y = COLLAPSE(X) is shorthand for S = COLLAPSE(X,1:ndims(X)).
%
%   Y = COLLAPSE(X,DIMS,FUN) accumulates the entries of T using the
%   accumulation function @FUN.
%
%   Examples
%   X = tenrand([4 4 4]);
%   Y = collapse(X,[2 3]) %<-- sum of entries in each mode-1 slice
%   Y = collapse(X,[1 2],@max) %<-- max entry in each mode-3 slice
%
%   See also TENSOR, TENSOR/SCALE.
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
% $Id: collapse.m,v 1.7 2010/03/19 23:46:30 tgkolda Exp $

if ~exist('fun', 'var')
    fun = @sum;
end
if ~exist('dims', 'var')
    dims = [];
end

dims = tt_dimscheck(dims,ndims(X));
remdims = setdiff(1:ndims(X),dims);

% Check for the case where we accumulate over *all* dimensions
if isempty(remdims)
    Y = fun(X.data(:));
    return;
end

% Calculate the size of the result
newsiz = size(X,remdims);

% Convert to a matrix where each row is going to be collapsed
A = double(tenmat(X,remdims,dims));

% Apply the collapse function
B = zeros(size(A,1),1);
for i = 1:size(A,1)
    B(i) = fun(A(i,:));
end

% Form and return the final result
Y = tensor(tenmat(B,1:numel(remdims),[],newsiz));




