function D = diff_matrix(n)
D = zeros(n,n);
d = zeros(1,n);
d(1) = -1;
d(2) = 1;
for j=1:n
    D(j,:)=circshift(d,[0,j-1]);
end
end


