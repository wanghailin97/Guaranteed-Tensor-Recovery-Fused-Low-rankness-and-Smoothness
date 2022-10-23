function [X,obj,err,iter] = LRMC_NN(M,omega,opts)

% Solve the Low-Rank Matrix Completion (LRMC) problem by ADMM
%
% min_X ||X||_*, s.t. P_Omega(X) = P_Omega(M)
%
% ---------------------------------------------
% Input:
%       M      -    d*n matrix
%       omega   -    index of the observed entries
%       lambda  -    >=0, parameter
%       opts    -    Structure value in Matlab. The fields are
%           opts.tol        -   termination tolerance
%           opts.max_iter   -   maximum number of iterations
%           opts.mu         -   stepsize for dual variable updating in ADMM
%           opts.max_mu     -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase mu
%           opts.detail      -   0 or 1
%
% Output:
%       X      -    d*n matrix
%       obj     -    objective function value
%       err     -    residual
%       iter    -    number of iterations
%
% version 1.0 - 22/06/2016
%
% Written by Canyi Lu (canyilu@gmail.com)
% 

tol = 1e-8; 
max_iter = 500;
rho = 1.1;
mu = 1e-4;
max_mu = 1e10;
detail = 1;

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'detail');       detail = opts.detail;          end

dim = size(M);
X = zeros(dim);
E = zeros(dim);
Lambda = zeros(dim);

%% main loop
iter = 0;
while iter<max_iter
    iter = iter + 1;  
    Xk = X;
    Ek = E;
    % update X
    [X,nuclearnormX] = prox_nuclear(-(E-M+Lambda/mu),1/mu);
    % update E
    E = -(X-M+Lambda/mu);
    E(omega) = 0;
    
    dY = X+E-M;  
    chgX = max(max(abs(Xk-X)));
    chgE = max(max(abs(Ek-E)));
    chg = max([chgX chgE max(abs(dY(:)))]);
    if detail
        if iter == 1 || mod(iter, 10) == 0
            obj = nuclearnormX;
            err = norm(dY,'fro');
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)]); 
        end
    end
    
    if chg < tol
        break;
    end 
    Lambda = Lambda + mu*dY;
    mu = min(rho*mu,max_mu);    
end
obj = nuclearnormX;
err = norm(dY,'fro');
