function X = my_ttm(X,V,dims,sizeOld,sizeNew,Ndim)
%X 只用是double的就好了
%V 要求是所有维度上都有，不乘的维度可以用[],如X = my_ttm(X,{U1,[],U3},(1,3),sizeX).
% sizeOld就是X的维数，其实很好算，sizeNew其实就是V里面的维度也很好算，不算是为了省时。

for n = dims(1:end) 
    order = [n,1:n-1,n+1:Ndim];
    temp = reshape(permute(X,order), sizeOld(n), []);
    temp = V{n}*temp;
    sizeOld(n) = sizeNew(n);    
    X = ipermute(reshape(temp, [sizeOld(n),sizeOld(1:n-1),sizeOld(n+1:Ndim)]),order);
end