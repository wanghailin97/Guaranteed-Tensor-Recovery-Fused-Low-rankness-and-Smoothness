function [L,obj,tsvd_rank] =HTNN_FFT(X,Omega,opts)

%Low-Rank Order-D(D>3) Tensor Completion  under  FFT 

% --------------------------------------------
% Input:
%                  X      -    Order-D tensor
%           Omega   -    index of the observed entries
%               opts    -    Structure value in Matlab. The fields are
%           opts.tol             -   termination tolerance
%           opts.max_iter   -   maximum number of iterations
%           opts.mu            -   stepsize for dual variable updating in ADMM
%           opts.max_mu   -   maximum stepsize
%           opts.rho           -   rho>=1, ratio used to increase mu
%           opts.DEBUG     -   0 or 1
%
% Output:
%       L                   -    Order-D tensor
%       obj               -    Objective function value
%      tsvd_rank     -    T_svd rank
%
% Written by  Wenjin Qin  (qinwenjin2021@163.com)
%

mu = 1e-3;
max_mu = 1e8;
tol = 1e-6; 
max_iter = 500;
rho = 1.2;
DEBUG = 1;

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'tol');         tol = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;    end
if isfield(opts, 'rho');         rho = opts.rho;              end
if isfield(opts, 'mu');          mu = opts.mu;                end
if isfield(opts, 'max_mu');      max_mu = opts.max_mu;        end
if isfield(opts, 'DEBUG');       DEBUG = opts.DEBUG;          end


dim = size(X);
p = length(size(X));
n = zeros(1,p);
for i = 1:p
    n(i) = size(X,i);
end

BarOmega = ones(dim) - Omega;
E = zeros(dim);
Y = E;
L = E;

for iter = 1 : max_iter
    Lk = L;
    Ek = E;
    
    % updata Mbar
    Mbar = L + E + Y/mu;
    Mbar = X.*Omega + Mbar.*BarOmega;
    
    % update L
    [L,tnnX,trank] =prox_htnn_F(Mbar-E-Y/mu,1/mu); 
    
    % updata M
    M = L + E + Y/mu;
    M = X.*Omega + M.*BarOmega;
    
    % update E
    E = (M-L-Y/mu).*BarOmega;
   
 
    dY = L+E-M;    
    chgX = max(abs(Lk(:)-L(:)));
    chgE = max(abs(Ek(:)-E(:)));
    chg = max([chgX chgE max(abs(dY(:)))]);
    if DEBUG
        if iter == 1 || mod(iter, 10) == 0
            obj = tnnX;
            err = chg;
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)]); 
        end
    end
    
    if chg < tol
        break;
    end 
    Y = Y + mu*dY;
    mu = min(rho*mu,max_mu);    
end
obj = tnnX;
tsvd_rank=trank;