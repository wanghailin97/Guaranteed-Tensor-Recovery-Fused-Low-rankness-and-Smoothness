%----- Correlted Total Viariation based Tensor Robust Principle Component Analysis -----%
function [X, E, obj, err, iter] = TCTV_TRPCA(M, opts)
% Solve the p-order Tensor Robust Principle Component Analysis via Tensor Correlted Total Viariation(TCTV) norm minimization by ADMM
% the transform in high-order TSVD uses DFT (default)
%
% min_{X \in R^{n1*n_2*...*n_d}} ||X||_tctv + lambda*||E||_1 s.t  M = X+ E
%
% ---------------------------------------------
% Input:
%       M       -    any p-order observed tensor
%       opts    -    Structure value in Matlab. The fields are
%           opts.directions          -   considered local smoothness along certain directions 
%           opts.transform           -   the transform case of TSVD, DFT, DCT and other invertible linear transform 
%           opts.transform_matrices  -   the transform matrices of TSVD for generalized invertible linear transform           
%           opts.tol                 -   termination tolerance
%           opts.max_iter            -   maximum number of iterations
%           opts.mu                  -   stepsize for dual variable updating in ADMM
%           opts.max_mu              -   maximum stepsize
%           opts.rho                 -   rho>=1, ratio that is used to increase mu
%           opts.detail             -   0 or 1, show the update details or not 
%
% Output:
%       X      -    recovered order-p tensor
%       E      -    corresponding sparse tensor
%       obj    -    objective function value
%       err    -    residual 
%       iter   -    number of iterations
%
% version 1.0 - 10/22/2022
%
% Written by Hailin Wang(wanghailin97@163.com) 
% 

%% default paremeters setting 
dim = size(M);
d   = ndims(M);

transform  = 'DFT';
for i = 3:d
transform_matrices{i-2} = dftmtx(dim(i)); 
end
lambda     = 1/sqrt(prod(dim)/min(dim(1),dim(2)));
directions = 1:d; 
tol        = 1e-8; 
max_iter   = 500;
rho        = 1.1;
mu         = 1e-4;
max_mu     = 1e10;
detail     = 1;

if ~exist('opts', 'var')
    opts = [];
end   
if isfield(opts, 'transform');          transform          = opts.transform;          end
if isfield(opts, 'transform_matrices'); transform_matrices = opts.transform_matrices; end
if isfield(opts, 'directions');         directions         = opts.directions;         end
if isfield(opts, 'tol');                tol                = opts.tol;                end
if isfield(opts, 'max_iter');           max_iter           = opts.max_iter;           end
if isfield(opts, 'rho');                rho                = opts.rho;                end
if isfield(opts, 'mu');                 mu                 = opts.mu;                 end
if isfield(opts, 'max_mu');             max_mu             = opts.max_mu;             end
if isfield(opts, 'detail');             detail             = opts.detail;             end

%% variables initialization
n = length(directions);
X        = randn(dim);
E        = zeros(dim);
Lambda   = zeros(dim);
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
    %% Update X -- solve TV by FFT 
	H = zeros(dim);
    for i = 1:n
       index = directions(i);
       H = H + porder_diff_T(mu*G{index}-Gamma{index},index); 
    end
    X = real( ifftn( fftn( mu*(M-E)+Lambda+H)./(mu*(1+T)) ) );
  
    %% Updata Gi -- proximal operator of TNN
    for i = 1:n
        index = directions(i);
        switch transform
            case 'DFT'
                [G{index},tnn_G{index}] = prox_htnn_F(porder_diff(X,index)+Gamma{index}/mu,1/(n*mu)); 
            case 'DCT'
                [G{index},tnn_G{index}] = prox_htnn_C(porder_diff(X,index)+Gamma{index}/mu,1/(n*mu));
            case 'other' % note: one need input the transform matrices
                [G{index},tnn_G{index}] = prox_htnn_C(transform_matrices,porder_diff(X,index)+Gamma{index}/mu,1/(n*mu));
        end
    end
    
    %% Update E 
	E = prox_l1(M-X+Lambda/mu, lambda/mu);
    
    %% Stop criterion
    dY   = M-X-E;    
    chgX = max(abs(Xk(:)-X(:)));
    chgE = max(abs(Ek(:)-E(:)));
    chg  = max([chgX chgE max(abs(dY(:)))]);
    if chg < tol
        break;
    end 
    
    %% Update detail display
    if detail
        if iter == 1 || mod(iter, 10) == 0
            obj = sum(cell2mat(tnn_G))/n;
            err = norm(dY(:),'fro');
            disp(['iter ' num2str(iter) ', mu=' num2str(mu) ...
                    ', obj=' num2str(obj) ', err=' num2str(err)]); 
        end
    end
    
    %% Update mulipliers: Lambda, Gamma, and mu
    Lambda = Lambda+mu*dY;
    for i = 1:n
        index = directions(i); 
        Gamma{index} = Gamma{index}+mu*(porder_diff(X,index)-G{index});
    end
    mu = min(rho*mu,max_mu);    
    
end
end

