function [P, P0, output] = cp_opt(Z,R,varargin)
%CP_OPT Fits a CP model to a tensor via optimization.
%
%   K = CP_OPT(X,R) fits an R-component CANDECOMP/PARAFAC (CP) model
%   to the tensor X. The result K is a ktensor. The function being
%   optimized is F(K) = 1/2 || X - K ||^2.
%
%   K = CP_OPT(X,R,'param',value,...) specifies additional
%   parameters for the method. Specifically...
%
%   'alg' - Specfies optimization algorithm (default: 'ncg')
%      'ncg'   Nonlinear Conjugate Gradient Method
%      'lbfgs' Limited-Memory BFGS Method
%      'tn'    Truncated Newton
%
%   'init' - Initialization for factor matrices. (default:
%   'random'). This can be a cell array with the initial matrices or
%   one of the following strings:
%      'random' Randomly generated via randn function
%      'nvecs'  Selected as leading left singular vectors of X(n)
%
%   'alg_options' - Parameter settings for selected optimization
%   algorithm. For example, type OPTIONS = NCG('defaults') to get
%   the NCG algorithm options which can then me modified as passed
%   through this function to NCG.
%
%   [K, U0] = CP_OPT(...) also returns the initial guess.
%
%   [K, U0, OUT] = CP_OPT(...) also returns a structure with the
%   optimization exit flag, the final relative fit, and the full
%   output from the optimization method.
%
%   See also CP_ALS, CP_FUN, TENSOR, SPTENSOR, KTENSOR.
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
% $Id: cp_opt.m,v 1.4 2010/03/22 18:07:17 tgkolda Exp $

%% Check for POBLANO
if ~exist('poblano_params','file')
    error(['CP_OPT requires Poblano Toolbox for Matlab. This can be ' ...
           'downloaded at http://software.sandia.gov/trac/poblano.']);
end

%% Error checking
if ~isa(Z,'tensor')
    Z = tensor(Z);
end

if (nargin < 2)
    error('Error: invalid input arguments');
end

%% Set parameters
params = inputParser;
params.addParamValue('alg', 'ncg', @(x) ismember(x,{'ncg','tn','lbfgs'}));
params.addParamValue('init', 'random', @(x) (iscell(x) || ismember(x,{'random','nvecs'})));
params.addOptional('alg_options', '', @isstruct);
params.parse(varargin{:});

%% Set up optimization algorithm
switch (params.Results.alg)
    case 'ncg'
        fhandle = @ncg;
    case 'tn'
        fhandle = @tn;
    case 'lbfgs'
        fhandle = @lbfgs;
end

%% Set up optimization algorithm options
if isempty(params.Results.alg_options)
    options = feval(fhandle, 'defaults');
else
    options = params.Results.alg_options;
end
        

%% Initialization
sz = size(Z);
N = length(sz);

if iscell(params.Results.init)
    P0 = params.Results.init;
elseif strcmpi(params.Results.init,'random')
    P0 = cell(N,1);
    for n=1:N
        P0{n} = randn(sz(n),R);
        for j=1:R
            P0{n}(:,j) = P0{n}(:,j) / norm(P0{n}(:,j));
        end
    end
elseif strcmpi(params.Results.init,'nvecs')
    P0 = cell(N,1);
    for n=1:N
        P0{n} = nvecs(Z,n,R);
    end
else
    error('Initialization type not supported')
end

%% Fit CP using CPOPT
normsqr = norm(Z)^2;
out = feval(fhandle, @(x)cp_fun(x,Z,normsqr), fac_to_vec(P0), options);

% compute factors and model fit
P = ktensor(cp_vec_to_fac(out.X, Z));
if nargout > 2
    output.ExitFlag  = out.ExitFlag;
    output.Fit = 100 * (1 - out.F /(0.5 * normsqr));
    output.OptOut = out;
end

%% Clean up final result
% Arrange the final tensor so that the columns are normalized.
P = arrange(P);
% Fix the signs
P = fixsigns(P);
