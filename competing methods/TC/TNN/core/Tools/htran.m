function B = htran(A,Mat_L)

%The conjugate transpose of a order-d(d>3) tensor
% Written by  Wenjin Qin  (qinwenjin2021@163.com)


p = length(size(A));
n = size(A); m= size(A);
m(2)=n(1);  m(1)=n(2);
L=1;
for i = 3:p
    L = L * n(i);
    m(i)=n(i);
end
A_Linear = lineartransform(A,Mat_L);
A_Tran=zeros(m);
for j = 1 : L
    A_Tran(:,:,j) = (A_Linear(:,:,j))';
end

B=inverselineartransform(A_Tran,Mat_L);

  