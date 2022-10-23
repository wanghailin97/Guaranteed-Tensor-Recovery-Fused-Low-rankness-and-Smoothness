function C = and(A,B)
%AND Logical AND (&) for sptensors.
%
%   See also SPTENSOR.
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
% $Id: and.m,v 1.6 2010/03/19 23:46:30 tgkolda Exp $

%% Observations for sparse matrix case.
% The result of a & 5 is sparse.
% The result of a & 0 is sparse.
% The result of a & full(a) is sparse.

%% Case 1: One argument is a scalar
if isscalar(B)
    if B == 0        
        C = sptensor([],[],size(A));
    else
        C = sptensor(A.subs,true,size(A));
    end   
    return;
end

% Call back with the arguments reversed.
if isscalar(A)
    C = and(B,A);
    return;
end

%% Case 2: Both x and y are tensors of some sort
% Check that the sizes match
if ~isequal(size(A),size(B))
    error('Must be tensors of the same size');
end

if isa(A,'sptensor') && isa(B,'sptensor')
    C = sptensor([A.subs; B.subs], [A.vals; B.vals], size(A), ...
        @(x) length(x) == 2);
    return;
end

if isa(B,'tensor')
    BB = sptensor(A.subs,B(A.subs),size(A));
    C = and(A,BB);
    return;    
end

%% Otherwise
error('The arguments must be two sptensors or an sptensor and a scalar.');
