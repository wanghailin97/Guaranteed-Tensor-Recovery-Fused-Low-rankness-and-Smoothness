clear all;
addpath(genpath(cd));
load Akiyo_RGB.mat
XT=Akiyo;
[n1 n2 n3 n4]=size(XT);


M{1}=sqrt(n3)*dctmtx(n3);
M{2}=sqrt(n4)*dctmtx(n4);
[U1,S1,V1]=ht_svd(XT,M);
s=diag(S1(:,:,1));
r=sum(s>0.01*max(s)); 
Xhat=zeros(n1,n2,n3,n4);
for i=1:r
   Xhat= Xhat + htprod_U( htprod_U(U1(:,i,:,:),S1(i,i,:,:),M),  htran(V1(:,i,:,:),M),M);
end
psnr_index=PSNR(XT, Xhat,max(XT(:)));
imshow( Xhat(:,:,:,50)/255)



