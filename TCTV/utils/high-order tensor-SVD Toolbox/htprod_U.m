function C = htprod_U(A,B,U)

% Tensor-tensor product of two order-d  tensors under generalized invertible linear transform

% A - Order-p tensor£ºn1*n2*n3*¡­¡­*np
% B - Order-p tensor: n2*l*n3*¡­¡­*np
% C - Order-p tensor: n1*l*n3*¡­¡­*np
%
% Written by  Wenjin Qin  (qinwenjin2021@163.com)

p = length(size(A));
n = zeros(1,p);
m = zeros(1,p);
for i = 1:p
    n(i) = size(A,i);
    m(i) = size(B,i);
end

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

sz = size(A); sz(2) = size(B,2);
C = zeros(sz);

L = ones(1,p);
for i = 3:p
    A=tmprod(A,U{i-2},i);
    B=tmprod(B,U{i-2},i);
    L(i) = L(i-1) * n(i);
end

for i = 1:L(p) 
    C(:,:,i)=A(:,:,i)*B(:,:,i);
end


for i = p:-1:3
    Mat_U=U{i-2};
    C = tmprod(C,inv(Mat_U),i);
end

end

