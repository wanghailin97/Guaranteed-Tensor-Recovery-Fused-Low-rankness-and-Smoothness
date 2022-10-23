function Approx = Approx_k(X,k)
dim = size(X);
[U, S, V] = tsvd(X, 'full');
Approx = zeros(dim);
for i = 1:k
    temp = tprod(tprod(U(:,i,:),S(i,i,:)),tran(V(:,i,:)));
    Approx = Approx + temp;
end
end







