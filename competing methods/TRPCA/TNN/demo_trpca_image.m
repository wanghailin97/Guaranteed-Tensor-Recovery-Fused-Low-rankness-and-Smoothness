addpath(genpath(cd))
clear

pic_name = './testimg.jpg';
X = double(imread(pic_name));

X = X/255;
maxP = max(abs(X(:)));
[n1,n2,n3] = size(X);
Xn = X;
rhos = 0.3;
ind = find(rand(n1*n2*n3,1)<rhos);
Xn(ind) = rand(length(ind),1);

opts.mu = 1e-4;
opts.tol = 1e-5;
opts.rho = 1.1;
opts.max_iter = 500;
opts.DEBUG = 1;

%% TNN-TRPCA
[n1,n2,n3] = size(Xn);
lambda = 1/sqrt(max(n1,n2)*n3);
[Xhat,E,err,iter] = trpca_tnn(Xn,lambda,opts);
Xhat = max(Xhat,0);
Xhat = min(Xhat,maxP);
psnr = PSNR(X,Xhat,maxP);

%% CTV-TRPCA
opts.detail = 1;
lambda = 1/sqrt(max(n1,n2)*n3);
[Xhat1,E,err,iter] = TRPCA_CTV(Xn,lambda,opts);
Xhat1 = max(Xhat1,0);
Xhat1 = min(Xhat1,maxP);
psnr1 = PSNR(X,Xhat1,maxP);


figure(1)
subplot(2,2,1)
imshow(X/max(X(:)))
subplot(2,2,2)
imshow(Xn/max(Xn(:)))
subplot(2,2,3)
imshow(Xhat/max(Xhat(:)))
subplot(2,2,4)
imshow(Xhat1/max(Xhat1(:)))




