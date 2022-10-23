function [U,S,V] = htsvd_fft(X)


% Order-d  tensors Singular Value Decomposition  under FFT
% Written by  Wenjin Qin  (qinwenjin2021@163.com)


 Ndim = length(size(X));
 Nway= zeros(1,Ndim);
for ii = 1:Ndim
    Nway(ii) = size(X,ii);
end

s1 = Nway; s1(2) = Nway(1);
s2 = Nway; s2(1) = Nway(2);

U=zeros(s1);S=zeros(Nway);V=zeros(s2);
L = ones(1,Ndim);
for ii_ = 3:Ndim
    X = fft(X,[],ii_);
    L(ii_) = L(ii_-1) * Nway(ii_);
end

[U(:,:,1),S(:,:,1),V(:,:,1)] = svd(X(:,:,1));

for j = 3 : Ndim
    for i = L(j-1)+1 : L(j)
      
        I = unfoldi(i,j,L);
        halfnj = floor(Nway(j)/2)+1;
        if I(j) <= halfnj && I(j) >= 2
              [U(:,:,i),S(:,:,i),V(:,:,i)] = svd(X(:,:,i));
        elseif I(j) > halfnj
            n_ = nc(I,j,Nway);
            i_ = foldi(n_,j,L);
            U(:,:,i) = conj(U(:,:,i_));
            V(:,:,i) = conj(V(:,:,i_));
            S(:,:,i) = conj(S(:,:,i_));
        end
    end
end


for jj = Ndim:-1:3
    U = ifft(U,[],jj);
    S = ifft(S,[],jj);
    V = ifft(V,[],jj);
end

end