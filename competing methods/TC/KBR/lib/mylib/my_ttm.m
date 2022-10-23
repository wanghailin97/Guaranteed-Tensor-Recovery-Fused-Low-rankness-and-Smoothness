function X = my_ttm(X,V,dims,sizeOld,sizeNew,Ndim)
%X ֻ����double�ľͺ���
%V Ҫ��������ά���϶��У����˵�ά�ȿ�����[],��X = my_ttm(X,{U1,[],U3},(1,3),sizeX).
% sizeOld����X��ά������ʵ�ܺ��㣬sizeNew��ʵ����V�����ά��Ҳ�ܺ��㣬������Ϊ��ʡʱ��

for n = dims(1:end) 
    order = [n,1:n-1,n+1:Ndim];
    temp = reshape(permute(X,order), sizeOld(n), []);
    temp = V{n}*temp;
    sizeOld(n) = sizeNew(n);    
    X = ipermute(reshape(temp, [sizeOld(n),sizeOld(1:n-1),sizeOld(n+1:Ndim)]),order);
end