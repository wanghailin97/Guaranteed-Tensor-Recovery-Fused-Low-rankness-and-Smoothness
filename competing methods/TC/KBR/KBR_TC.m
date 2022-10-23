function [ X, C, U ] = KBR_TC(T, Omega, Par)
% sloving following KBR-based tensor completion problem.
% min_X: P_ls( C ) + lambda Prod_j{(P_ls^*( M_j(j))} 
% s.t.:  X_Omega = T_Omega, ttm(C,{U_1,U_2,...,U_N}) -X = 0, X-M_j = 0, forall j =  1,2,..,N;
%
% Input arguments:
%   T      ...... the corrupted tensor. Please make sure T is in range [0, 1].
%   Omega  ...... a mask of the observed elements of T
%   Par    ...... an option structure whose fields are as follows:
%      lambda ... the compromise parameter in ITS, usually setted in [0.1,10];
%      mu     ... initial mu in ADMM algorithm;
%      rho    ... parameter control the increasing speed of mu
%      maxIter .. max iteration step number
%      tol     .. termination threshold
%
% Output arguments:
%   X     ......  the reconstruct tensor
%   C     ......  core tensor obtained by ITS-TC
%   U     ......  Bases obtained by ITS-TC
%
% more detail can be found in [1]
% [1] Qi Xie, Qian Zhao, Deyu Meng*, Zongben Xu. 
%     Kronecker-Basis-Representation Based Tensor Sparsity and Its Applications to Tensor Recovery,  
%     IEEE Transactions on Pattern Analysis and Machine Intelligence, 2017.
%
%       @article{Qi2017Kronecker,  
%       title={Kronecker-Basis-Representation Based Tensor Sparsity and Its Applications to Tensor Recovery},  
%       author={Qi, Xie and Qian, Zhao and Meng, Deyu and Xu, Zongben},  
%       journal={IEEE Transactions on Pattern Analysis & Machine Intelligence},  
%       volume={PP},  
%       number={99},  
%       pages={1-1},  
%       year={2017},
%       }
%
% by Qi Xie 
% 10/24/2017
%==========================================================================

X         = T.*Omega;
X(~Omega) = mean(T(Omega));
sizeT     = size(T);
ndim      = length(sizeT);
dim1Xdim2 = circshift(sizeT, [1,1]).*circshift(sizeT, [2,2]);
if nargin<3
    MaxIter        = 100;
    tol            = 1e-7;
    rank           = min(sizeT,dim1Xdim2);
    mu             = 1;
    rho            = 1.05; 
    lambda         = 1;
else
    if isfield(Par,'maxIter')   ; MaxIter = Par.maxIter;       else MaxIter = 100;  end
    if isfield(Par,'tol')       ; tol = Par.tol;               else tol = 1e-7;     end
    if isfield(Par,'rank')      ; rank = Par.rank;             else rank = size(T); end
    if isfield(Par,'mu')        ; mu = Par.mu;                 else mu = 1;         end
    if isfield(Par,'rho')       ; rho = Par.rho;               else rho = 1.05;     end
    if isfield(Par,'lambda')    ; lambda = Par.lambda;         else lambda = 1;     end
    
end

normT     = norm(T(:));
Bet       = cell(ndim, 1);
M         = cell(ndim, 1);
tempX     = cell(ndim, 1);
sValue    = cell(ndim, 1);
Mnorm     = zeros(ndim, 1);
for i = 1:ndim
    M{i}      = X;
    Bet{i}    = zeros(sizeT);
    tempX{i}  = Unfold(X, sizeT, i);
    sValue{i} = svd(tempX{i}, 0);
    Mnorm(i)  = sum(sValue{i}>eps);
end
alpha     = circshift(Mnorm, [1,1]).*circshift(Mnorm, [2,2]); %computing the  weigth
sumAlpha  = sum(alpha);

mu        = mu*sumAlpha;
muLam     = mu/100;
muBet     = mu;
muGam     = mu;
Lam       = zeros(sizeT);
Gam       = zeros(sizeT);
Msum      = zeros(sizeT);
BetSum    = zeros(sizeT);
temp_n    = zeros(ndim,1);

if nargin<4
    [C, U]    = tensorSVD(X+Lam/muLam,rank);
    [C, ~]    = ClosedWL1(C,1/muLam,eps);
end

Y         = double(ttm(tensor(C),U));
%%main loop
% fprintf('iter         ')
for i = 1:MaxIter
    preX = X;
 
    % update U
    for k = 1:ndim
        unfoTemp    = Unfold(X+Lam/muLam, sizeT, k);
        tempC       = double(ttm(tensor(C),U,-k));
        UnfoldC     = Unfold( tempC, rank, k);
        A           = unfoTemp*UnfoldC';
        dim         = size(A);
        [V1,~,V2]   = svd(A);
        U{k}        = V1(:,1:dim(2))*V2';
        Ut{k}       = U{k}';
    end
    
    % update C
    C       = double(ttm(tensor(X+Lam/muLam),Ut));
    [C, n]  = ClosedWL1(C,1/muLam,eps);
    %         C = max(min(C+1/muLam,0),C-1/muLam);
    
    % update Y
    Y       = double(ttm(tensor(C),U));
    
    % update M
    Msum   = 0*Msum;
    BetSum = 0*BetSum;
    for k = 1:ndim
        [tempX{k}, temp_n(k), sValue{k}, Mnorm(k)] = Pro2logSum(Unfold(X + Bet{k}/muBet, sizeT, k), lambda*alpha(k)/muBet);
        M{k}      = Fold(tempX{k}, sizeT, k);
        Msum      = Msum + M{k};
        BetSum    = BetSum + Bet{k};
        alpha     = circshift(Mnorm, [1,1]).*circshift(Mnorm, [2,2]); %computing the  weigth
    end
    
    % update X
    X(~Omega)   = (muLam*Y(~Omega) - Lam(~Omega) +muBet*Msum(~Omega) - BetSum(~Omega))/(ndim*muBet+muLam);
    X(Omega)    = (muLam*Y(Omega) - Lam(Omega) + muBet*Msum(Omega) - BetSum(Omega)+ muGam*T(Omega) - Gam(Omega))/(ndim*muBet+muLam+muGam);
   
    % update multiplayer
    Lam        = Lam + muLam * (X - Y);
    for k = 1:ndim
        Bet{k} = Bet{k} + muBet*( X - M{k});
    end
    Gam(Omega) = Gam(Omega)+muGam*(X(Omega)-T(Omega));
    
    % update mu
    muLam = rho * muLam;
    muBet = rho * muBet;
    muGam = rho * muGam;
    
    %%stop conditon
    %     fprintf('\b\b\b\b\b%5i',i);
    if norm(preX(:)-X(:))/normT <tol
        break
    end
end

if i == MaxIter
    disp(['Does not convergence in  ' num2str(MaxIter) ' steps.'])
end

% fprintf('\n')
f= sum(abs(C(:))./(eps+abs(C(:))));
g= prod(Mnorm);
end


function [X, n, SigmaNew,wNorm] = Pro2logSum(Z, tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min: 1/2*||Z-X||^2 + (-log(eps)) * tau * P_ls^*(X)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% new
% [m, n] = size(Z);
% if m < n
    AAT = Z*Z';
    [S, Sigma, ~] = svd(AAT);
    Sigma         = sqrt(diag(Sigma));    
    tol           = eps;%max(size(Z)) * eps(max(Sigma));
%     [SigmaNew,n]  = ClosedWNNM(Sigma,tau,tol);
    temp      = (Sigma-tol).^2-4*(tau-tol*Sigma);
    ind       = find (temp>0);
    n         = length(ind);
    SigmaNew  = (Sigma(1:n)-tol+sqrt(temp(1:n)))/2;
    wNorm         = sum((log(SigmaNew+eps)+36.04)/(36.04+7)); % log(eps) = 36.0437
    SigmaNew      = SigmaNew ./ Sigma(1:n);
    X = S(:, 1:n) * diag(SigmaNew) * S(:, 1:n)' * Z;
%     return;
% else
%     [X, n, SigmaNew, wNorm] = Pro2logSum(Z', tau);
%     X = X';
%     return;
end
