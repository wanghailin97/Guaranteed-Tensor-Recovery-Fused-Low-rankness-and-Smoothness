function [ B, T] = KBR_RPCA(D,beta,gamma, Par, init)
% sloving following KBR-based tensor RPCA problem.
% min_X: |Cx|_w.1+lambda prod(|M_i|_w.*) +beta |T|_1 + gamma/2 |ttm(C,U)_(i)+T-X| _F^
% s.t.  ttm(C,U)_(i)_(i) - M_i =0;   
%
% Input arguments:
%   D      ...... the corrupted tensor. Please make sure T is in range [0, 1].
%   beta   ...... tuning parameter compromising the recovery tensor and sparse noise terms
%   gamma  ...... tuning parameter compromising the recovery tensor and L1 noise terms
%   Par    ...... an option structure whose fields are as follows:
%      lambda  .. the compromise parameter in ITS, usually setted in [0.1,10];
%      mu      .. initial mu in ADMM algorithm;
%      rho     .. parameter control the increasing speed of mu
%      maxIter .. max iteration step number
%      tol     .. termination threshold
%      Rank0   .. a up bound of rank, which could be set as the size of input tensor 
%   init    ......input initialization whose fields are as follows:
%      B       .. the initial high-order sparse component
%      T       .. the initial 1-order sparse component
%
% Output arguments:
%   B     ......  output high-order sparse component
%   T     ......  output 1-order sparse component
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

%% initialization

sizeD     = size(D);
ndim      = length(sizeD);
dim1Xdim2 = circshift(sizeD, [1,1]).*circshift(sizeD, [2,2]);
normD     = norm(D(:));
if nargin <2
    beta  = 2*sqrt(max(sizeD));
    gamma = beta*100;
elseif nargin <3
    gamma= beta*100;
end
if nargin<4
    MaxIter        = 1000;
    tol            = 1e-5;
    Rank0          = min(sizeD,dim1Xdim2);
    mu             = 10;
    rho            = 1.05; 
    lambda         = 0.1;
else
    if isfield(Par,'maxIter')   ; MaxIter = Par.maxIter;       else MaxIter = 1000;               end;
    if isfield(Par,'tol')       ; tol = Par.tol;               else tol = 1e-5;                   end;
    if isfield(Par,'Rank0')     ; Rank0 = Par.Rank0;           else Rank0 = min(sizeD,dim1Xdim2); end;
    if isfield(Par,'mu')        ; mu = Par.mu;                 else mu = 10;                      end;
    if isfield(Par,'rho')       ; rho = Par.rho;               else rho = 1.05;                   end;
    if isfield(Par,'lambda')    ; lambda = Par.lambda;         else lambda = 0.1;                 end;
end
if nargin<5
    initT      = zeros(sizeD);
    initB      = D;
else
    if isfield(init,'B')     ;     initB = init.B;       else     initB  = ones(sizeD)*mean(y);  end
    if isfield(init,'T')     ;     initT = init.T;       else     initT  = zeros(sizeD);     end
end

Rank           = round([0.5,0.5,0.1].*Rank0); % start from little rank for computational cost saving
%% initialization about M
M         = cell(ndim, 1);
LamM      = cell(ndim, 1);
tempB     = cell(ndim, 1);
sValue    = cell(ndim, 1);
Mnorm     = zeros(ndim, 1);
Msum      = zeros(sizeD);
Anum      = [2,2,8];
for i = 1:ndim
    M{i}      = initB;
    LamM{i}   = zeros(sizeD);
    tempB{i}  = Unfold(initB, sizeD, i);
    sValue{i} = svd(tempB{i}, 'econ');
    Mnorm(i)  = min(sum(sValue{i}>Anum(i)),Rank(i));
    Msum      = Msum + Fold(M{i},sizeD,i);
end
weigh     = circshift(Mnorm, [1,1]).*circshift(Mnorm, [2,2]); %computing the  weigth
mu        = mu*(1+lambda);
beta      = beta*(1+lambda);
gamma     = gamma*(1+lambda);
%% initialization about C
[C,U]   = tensorSVD2(initB,Rank);
Ut      = cell(ndim,1);
[C, ~]  = ClosedWL1(C,2/beta,eps);
% B       = my_ttm(C,U,1:ndim,Rank,sizeD,ndim);

%% initialization about other parameters
B         = initB;
T         = initT;
% F         = zeros(sizeD,ndim);
LamMsum   = zeros(sizeD);
temp_n    = zeros(1,ndim);
%% main loop
fprintf('iter£º         ')
for i = 1:MaxIter
    
    preB = B;
    
    %% update T 
    T = min(max(D-B-beta/gamma,0),D-B+beta/gamma);
    
    %% update U C
    for j = 1:ndim
        unfoTemp    = Unfold((gamma*(D-T)+mu*Msum-LamMsum)/(gamma+ndim*mu), sizeD, j);
        tempC       = my_ttm(C,U,[1:j-1,j+1:ndim],Rank,sizeD,ndim);
        UnfoldC     = Unfold( tempC, Rank, j);
        tempMatix   = unfoTemp*UnfoldC';
        [V1,~,V2]   = svd(tempMatix,'econ');
        U{j}        = V1*V2';
        Ut{j}       = U{j}';
    end
    C       = my_ttm((gamma*(D-T)+mu*Msum-LamMsum)/(gamma+ndim*mu),Ut,1:ndim,sizeD,Rank,ndim);
    [C, ~]  = ClosedWL1(C,1/mu,eps);
    %% computing B
    B       = my_ttm(C,U,1:ndim,Rank,sizeD,ndim);
    
    %% update M
    Msum    = 0*Msum;
    LamMsum = 0*LamMsum;
    for j = ndim:-1:1
        [tempB{j}, temp_n(j), sValue{j}, Mnorm(j)] = Pro2logSum(Unfold(B + LamM{j}/mu, sizeD, j), lambda*weigh(j)/mu);%, Rank0(j));
        M{j}      = Fold(tempB{j}, sizeD, j);
        Msum      = Msum + M{j};
        weigh     = circshift(Mnorm, [1,1]).*circshift(Mnorm, [2,2]); %computing the  weigth
        LamM{j}   = LamM{j}+mu*(B-M{j}); % update multiplayer
        LamMsum   = LamMsum + LamM{j};  % update multiplayer
    end  
 
    %% update mu
    mu        = mu*rho;
    
    [C,U,Rank] = ChangerRank(C,U,Rank,temp_n,Rank0,ndim,sizeD); % changing rank, just for computational cost saving
    
    stopCond = norm(preB(:)-B(:))/normD ;

    %% stop conditon
    fprintf('\b\b\b\b\b%5i',i);
% if mod(i,10)==1
%     figure(1);
%     subplot(221);imshow(B(:,:,1));title('B');
%     subplot(222);imshow(T(:,:,1),[]);
%     subplot(223);imshow(M{3}(:,:,1));title('M');
%     subplot(224);imshow(D(:,:,1)-B(:,:,1)-T(:,:,1),[]);
%     disp(['iter = ' num2str(i) ';    tol = ' num2str(stopCond) '.'])
%     disp(['      rankM = ' num2str(temp_n) ';    rankC = ' num2str(Rank) '.'])
% end
    if stopCond <tol
        break
    end
end
fprintf('\n')
if i == MaxIter
    disp(['Does not convergence in  ' num2str(MaxIter) ' steps.'])
end
f= sum(abs(C(:))./(eps+abs(C(:))));
g= prod(Mnorm);
end

function [C,U,Rank] = ChangerRank(C,U,Rank,tempR,Rank0,ndim,sizeD)
IndChange  = find(tempR>Rank);
newR       = min(Rank+1,tempR+20);
for i = IndChange(1:end)
    newR(i) = round(min(Rank(i)+0.05*Rank0(i),Rank0(i)));
end
if ~isequal(newR,Rank)
    if sum(newR>Rank)>0
        IndBiger  = find(newR>Rank);
        for i = IndBiger(1:end)
            temp = zeros(sizeD(i),newR(i));
            temp(:,1:Rank(i)) = U{i};
            U{i} = temp;
        end
        ind       = true(Rank);
        Rank      = max(newR,Rank);
        temp      = zeros(Rank);
        temp(ind) = C;
        C         = temp;
    end
    %% sort the element of a tensor
    for i = 1:ndim
        tempC      = Unfold(C,Rank,i);
        tempNorm   = sum(tempC.^2,2);
        [~,Ind] = sort(tempNorm,'descend');
        tempC      = tempC(Ind(1:newR(i)),:);
        Rank(i)    = newR(i);
        C          = Fold(tempC,Rank,i);
        U{i} = U{i}(:,Ind(1:newR(i)));
    end
end
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
