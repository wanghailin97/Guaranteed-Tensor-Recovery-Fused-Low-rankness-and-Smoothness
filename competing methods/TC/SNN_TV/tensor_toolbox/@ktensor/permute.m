function b = permute(a,order)
%PERMUTE Permute dimensions of a ktensor.
%
%   B = PERMUTE(A,ORDER) rearranges the dimensions of A so that they
%   are in the order specified by the vector ORDER. The tensor
%   produced has the same values of A but the order of the subscripts
%   needed to access any particular element are rearranged as
%   specified by ORDER.
%
%   See also KTENSOR.
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
% $Id: permute.m,v 1.8 2010/03/19 23:46:30 tgkolda Exp $

N = ndims(a);

if ~isequal(1:N,sort(order))
  error('Invalid permuation');
end

b = ktensor(a.lambda, a.u(order));




