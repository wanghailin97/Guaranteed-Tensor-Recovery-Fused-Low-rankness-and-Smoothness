function Approx = Approx_k_U(ModeX,X,k)
dim = size(X);
[U, S, V] = htsvd(X,ModeX);
Approx = zeros(dim);
for i = 1:k
    temp = htprod_U(htprod_U(U(:,i,:),S(i,i,:),ModeX),htran(V(:,i,:),ModeX),ModeX);
    Approx = Approx + temp;
end
end







