function X = tenzeros(sz,varargin)
%TENZEROS Zeros tensor.
%
%   X = TENZEROS(SZ) forms a tensor of size SZ with all zeros.
%
%   X = TENZEROS(ORD,DIM) forms a tensor of order ORD with each mode being
%   of size D.
%
%   TENZEROS(SZ) is equivalent to TENSOR(ZEROS(SZ(1),SZ(2),...),SZ).
%
%   See also TENSOR, ZEROS.
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
% $Id: tenzeros.m,v 1.7 2010/03/19 23:46:31 tgkolda Exp $

if isempty(sz)
    X = tensor;
    return;
end

if nargin == 2
    order = sz;
    dim = varargin{1};
    sz = dim * ones(1,order);
end

data = zeros([sz 1 1]);
X = tensor(data,sz);

