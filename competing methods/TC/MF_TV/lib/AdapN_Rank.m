%%%compute the n-rank
function rank=AdapN_Rank(X,ratio)
rank=zeros(1,3);
for i = 1:3
a1 = tenmat(X,i);
a1 = double(a1);
temp = a1*a1';
s = eig(temp);
s=abs(sqrt(s));
big=s(end);
rank(i)=size(find(s>=ratio*big),1);
end
