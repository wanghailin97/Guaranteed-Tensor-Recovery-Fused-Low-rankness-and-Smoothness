function Z = plus(X,Y)
%PLUS Binary addition (+) for tenmat. 
%
%   See also TENMAT, TENMAT/TENMATFUN.
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
% $Id: plus.m,v 1.7 2010/03/19 23:46:30 tgkolda Exp $

fun = @plus;

% One argument is a scalar
if ((prod(size(X)) == 1 || prod(size(Y)) == 1))
  if (prod(size(Y)) == 1) && isa(X,'tenmat')
    Z = X;
    Z.data = fun(Z.data,Y);
  else
    Z = Y;
    Z.data = fun(X,Z.data);
  end
  return;
end


% Both arguments are tenmats
Z = tenmat(Y);
if ~(isequal(size(Y),size(Z)))
  error('Tenmat size mismatch.')
end
Z.data = fun(X.data,Z.data);
