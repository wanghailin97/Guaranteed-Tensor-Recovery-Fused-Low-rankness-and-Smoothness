function [X,tnn,trank] = prox_dcttnn(Y,rho)

% The proximal operator of the tensor nuclear norm of a 3 way tensor
%
% min_X rho*||X||_*+0.5*||X-Y||_F^2
%
% Y     -    n1*n2*n3 tensor
%
% X     -    n1*n2*n3 tensor
% tnn   -    tensor nuclear norm of X
% trank -    tensor tubal rank of X
%

% 

[n1,n2,n3] = size(Y);
X = zeros(n1,n2,n3);
Y = dct(Y,[],3);
tnn = 0;
trank = 0;
% first frontal slice
% [U,S,V] = svd(Y(:,:,1),'econ');
% S = diag(S);
% S = max(S-rho,0);
% r = length(find(S~=0));
% S = S(1:r);
% X(:,:,1) = U(:,1:r)*diag(S)*V(:,1:r)';
% tnn = tnn+sum(S);
% trank = max(trank,r);

% i=2,...,halfn3
% halfn3 = round(n3/2);
for i = 1 : n3
    [U,S,V] = svd(Y(:,:,i),'econ');
    S = diag(S);
    S = max(S-rho,0);
    r = length(find(S~=0));
    S = S(1:r);
    X(:,:,i) = U(:,1:r)*diag(S)*V(:,1:r)';    
%     X(:,:,n3+2-i) = conj(X(:,:,i));
    tnn = tnn+sum(S)*2;
    trank = max(trank,r);
% end
% 
% % if n3 is even
% if mod(n3,2) == 0
%     i = halfn3+1;
%     [U,S,V] = svd(Y(:,:,i),'econ');
%     S = diag(S);
%     S = max(S-rho,0);
%     r = length(find(S~=0));
%     S = S(1:r);
%     X(:,:,i) = U(:,1:r)*diag(S)*V(:,1:r)';
%     tnn = tnn+sum(S);
%     trank = max(trank,r);
end
tnn = tnn/n3;
X = idct(X,[],3);
