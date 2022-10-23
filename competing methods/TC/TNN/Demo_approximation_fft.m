clear all;
addpath(genpath(cd));
load Akiyo_RGB.mat
XT=Akiyo;
[n1 n2 n3 n4]=size(XT);



[U1,S1,V1]=htsvd_fft(XT);
s=diag(S1(:,:,1));
r=sum(s>0.01*max(s)); 
Xhat=zeros(n1,n2,n3,n4);
for i=1:r
   Xhat= Xhat + htprod_fft( htprod_fft(U1(:,i,:,:),S1(i,i,:,:)),  htran(V1(:,i,:,:),'fft'));
end
psnr_index=PSNR(XT, Xhat,max(XT(:)));
imshow( Xhat(:,:,:,50)/255)



