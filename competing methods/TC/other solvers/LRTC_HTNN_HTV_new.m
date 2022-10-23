%-----Low-tubal-rank Tensor Completion via TNN plus 3D Total Viariation-----%
function [X,obj,err,iter] =LRTC_HTNN_HTV_new(M,Omega,lambda,opts) 
% Solve the p-order Tensor Completion via Tensor Nuclear Norm(TNN) plus High-order Total Variation(HTV) norm minimization by ADMM
%
% min_{X \in R^{n1*n_2*...*n_d}} ||X||_*+\lambda*||X||_htv s.t. P_Omega(X) = P_Omega(M)
%
% ---------------------------------------------
% Input:
%       M       -    any p-order observed tensor
%       opts    -    Structure value in Matlab. The fields are
%           opts.directions -   consider nonlocal smoothness along certain directions
%           opts.tol        -   termination tolerance
%           opts.max_iter   -   maximum number of iterations
%           opts.mu         -   stepsize for dual variable updating in ADMM
%           opts.max_mu     -   maximum stepsize
%           opts.rho        -   rho>=1, ratio used to increase mu
%           opts.detail     -   0 or 1, show some update details of algorithm
%
% Output:
%       X       -    recovered order-p tensor
%       obj     -    objective function value
%       err     -    residual 
%       iter    -    number of iterations
%
% version 1.0 - 12/12/2021
%
% Written by Hailin Wang(wanghailin97@163.com) 
% 
%% paremeters setting
dim = size(M);
d   = ndims(M);

directions = 1:d; 
tol       = 1e-8; 
max_iter  = 500;
rho       = 1.1;
mu        = 1e-4;
max_mu    = 1e10;
detail    = 0;

if ~exist('opts', 'var')
    opts = [];
end    
if isfield(opts, 'directions');  directions = opts.directions;     end
if isfield(opts, 'tol');         tol      = opts.tol;              end
if isfield(opts, 'max_iter');    max_iter = opts.max_iter;         end
if isfield(opts, 'rho');         rho      = opts.rho;              end
if isfield(opts, 'mu');          mu       = opts.mu;               end
if isfield(opts, 'max_mu');      max_mu   = opts.max_mu;           end
if isfield(opts, 'detail');      detail   = opts.detail;           end

%% variables initialization
n = length(directions);
X        = randn(dim);
X(Omega) = M(Omega);
Y        = zeros(dim);
E        = zeros(dim);
Lambda1  = zeros(dim);
Lambda2  = zeros(dim);
for i = 1:n
    index        = directions(i);
    G{index}     = porder_diff(X,index); 
    Gamma{index} = zeros(dim); 
end

%% FFT setting
T = zeros(dim);
for i = 1:n
    Eny = diff_element(dim,directions(i));
    T   = T + Eny; 
end 

%% main loop
iter = 0;
while iter<max_iter
    iter = iter + 1;  
    Xk = X;
    Ek = E;
    Yk = Y;
    %% Update X -- solve TV by FFT 
    H = zeros(dim);
    for i = 1:n
       index = directions(i);
       H = H + porder_diff_T(mu*G{index}-Gamma{index},index); 
    end
    X = real( ifftn( fftn( mu*(M+Y-E)+Lambda1+Lambda2+H)./(mu*(2+T)) ) );
  
    %% Update Y -- proximal operator of TNN
    [Y,tnn_Y] = porder_prox_tnn(X-Lambda1/mu,1/mu); 
    
    %% Updata Gi -- soft-thresholding operator
    for i = 1:n
        index = directions(i);
        G{index} = prox_l1(porder_diff(X,index)+Gamma{index}/mu,lambda/mu); 
    end
    
    %% Update E 
    E          = M-X+Lambda2/mu;
    E(Omega)   = 0;
    
    %% Stop criterion
    dY   = M-X-E;    
    chgX = max(abs(Xk(:)-X(:)));
    chgE = max(abs(Ek(:)-E(:)));
    chgY = max(abs(Yk(:)-Y(:)));
    chg  = max([chgX chgE chgY max(abs(dY(:)))]);
    if chg < tol
        break;
    end 
    
    %% Update detail display
    if detail
        if iter == 1 || mod(iter, 10) == 0
            mat_G = cell2mat(G);
            obj = tnn_Y+lambda*norm(mat_G(:),1);
            err = norm(dY(:),'fro');
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)]); 
        end
    end
    
    %% Update mulipliers: Lambda, Gamma, and mu
    Lambda1 = Lambda1+mu*(Y-X);
    Lambda2 = Lambda2+mu*dY;
    for i = 1:n
        index = directions(i); 
        Gamma{index} = Gamma{index}+mu*(porder_diff(X,index)-G{index});
    end
    mu = min(rho*mu,max_mu);    
    
end
end

