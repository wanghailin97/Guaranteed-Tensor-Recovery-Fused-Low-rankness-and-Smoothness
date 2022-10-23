function Approx = Approx_k_4D(X,k)
dim = size(X);
[U, S, V] = htsvd_fft(X);
Approx = zeros(dim);
for i = 1:k
    temp = htprod_fft(htprod_fft(U(:,i,:,:),S(i,i,:,:)),htran(V(:,i,:,:),'fft'));
    Approx = Approx + temp;
end
end







