function [X,htnn,tsvd_rank] = prox_htnn_U(M,Y,rho)

% The proximal operator for the order-D tensor nuclear norm under generalized invertible linear transform
%
% you need input the specific invertible linear transforms M which stores the 3-th to p-th transform matrix in M{1} to M{p-2}
%
% Written by  Wenjin Qin  (qinwenjin2021@163.com)
%


p = length(size(Y));
n = zeros(1,p);
for i = 1:p
    n(i) = size(Y,i);
end
X = zeros(n);

L = ones(1,p);
for i = 3:p
     Y = tmprod(Y,M{i-2},i);
    L(i) = L(i-1) * n(i);
end

htnn = 0;
tsvd_rank = 0;
       
for i=1:L(p)
[U,S,V] = svd(Y(:,:,i),'econ');
S = diag(S);
r = length(find(S>rho));
if r>=1
    S =max( S(1:r)-rho,0);
    X(:,:,i) = U(:,1:r)*diag(S)*V(:,1:r)';
    htnn = htnn+sum(S);
    tsvd_rank = max(tsvd_rank,r);
end
end

rho=1;
for j=3:p
     Tran_M=M{j-2};
     a=sum(diag(Tran_M*(Tran_M)'))/n(j);
     rho=rho*a;
end

htnn = htnn/rho;

for i = p:-1:3
    X = tmprod(X,inv(M{i-2}),i);
end  

X = real(X);


