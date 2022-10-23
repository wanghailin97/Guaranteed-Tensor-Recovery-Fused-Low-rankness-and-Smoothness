function C = htprod_fft(A,B)

% Tensor-tensor product of two order-d  tensors under FFT

% A - Order-p tensor£ºn1*n2*n3*¡­¡­*np
% B - Order-p tensor: n2*l*n3*¡­¡­*np
% C - Order-p tensor: n1*l*n3*¡­¡­*np
%
% Written by  Wenjin Qin  (qinwenjin2021@163.com)



p = length(size(A));
n=size(A);
m=size(B);
for i = 3:p
    if n(i) ~= m(i)
        t = 1;
        break;
    else t = 0;
    end
end

if n(2) ~= m(1) || t == 1
    error('Inner tensor dimensions must agree.');
end

C = zeros([ size(A,1)  size(B,2)  n(3:end) ]);
L = ones(1,p);
for i = 3:p
    A = fft(A,[],i);
    B = fft(B,[],i);
    L(i) = L(i-1) * n(i);
end

C(:,:,1) = A(:,:,1)*B(:,:,1);   

for j = 3 : p
    for i = L(j-1)+1 : L(j)
        I = unfoldi(i,j,L);
        halfnj = floor(n(j)/2)+1;
        if I(j) <= halfnj && I(j) >= 2
            C(:,:,i) = A(:,:,i)*B(:,:,i);   
        elseif I(j) > halfnj
            n_ = nc(I,j,n);
            i_ = foldi(n_,j,L);
            C(:,:,i) = conj(C(:,:,i_));
        end
    end
end

for i = p:-1:3
    C = (ifft(C,[],i));
end
C = real(C);
end
