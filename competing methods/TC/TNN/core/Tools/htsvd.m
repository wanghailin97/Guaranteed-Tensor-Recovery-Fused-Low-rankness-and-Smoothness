function [U,S,V] = htsvd(X,UU)


% Order-d  tensors Singular Value Decomposition  under generalized invertible linear transform
% Written by  Wenjin Qin  (qinwenjin2021@163.com)

 Ndim = length(size(X));
 Nway= size(X);
 
s1 = Nway; s1(2) = Nway(1);
s2 = Nway; s2(1) = Nway(2);

U=zeros(s1);S=zeros(Nway);V=zeros(s2);
L = ones(1,Ndim);
for i = 3:Ndim
    % X = nmodeproduct(X,UU{i-2},i);
    X = tmprod(X,UU{i-2},i);
    L(i) = L(i-1) * Nway(i);
end

for i = 1 : L(Ndim)
   [U(:,:,i),S(:,:,i),V(:,:,i)] = svd(X(:,:,i));
end


for j = Ndim:-1:3
    U = tmprod(U,inv(UU{j-2}),j);
    S = tmprod(S,inv(UU{j-2}),j);
    V = tmprod(V,inv(UU{j-2}),j);
end

end