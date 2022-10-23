%%
clear all
addpath(genpath(cd));
randn('seed',2013);
randn('seed',2013);

%%  Face data
load('YaleFace.mat');
X = YaleFace./max(YaleFace(:));
[d1, d2 ,d3,d4] = size(X);
maxP = max(abs(X(:)));
maxP1 = max(YaleFace(:));
Nways=size(X);
     
%%  Set  Omega
sr = 0.1;  %% ssampling ratio \rho
fprintf('Sampling ratio = %0.8e\n',sr);
temp = randperm(prod(Nways));
kks = round((sr)*prod(Nways));
mark = zeros(Nways); 
mark(temp(1:kks)) = 1;
Y=X.*mark;

shift_dim=[3,4,1,2];
XT = permute(X,shift_dim);
mark = permute(mark,shift_dim);
[n1, n2 ,n3,n4] = size(XT);

%% initial parameters
opts.mu = 1e-4;
opts.max_mu = 1e8;
opts.max_iter =500;
opts.DEBUG = 1;
opts.rho = 1.2;
opts.tol = 1e-6;


%% Fast  Fouire   Transform (FFT)
  fprintf('===== t-SVD by FFT =====\n');
 
   t0=tic;
     [Xhat,~,~ ] =HTNN_FFT(XT,mark,opts);
   time = toc(t0);
   Xhat1=max(0,Xhat);
   Xhat2=min(maxP,Xhat1);
   Xhat2=permute(Xhat2,shift_dim);

%% print the relative error, psnr
    Error = norm(Xhat2(:)-X(:))/norm(X(:));
    fprintf('Relative error = %0.8e\n',Error);
    psnr_index = PSNR(Xhat2,X,maxP);
    [~,ssim,fsim]=quality(Xhat2*maxP1,X*maxP1);
    fprintf('PSNR = %0.8e\n',psnr_index);



    figure(1);
    imshow(X(:,:,11,10))
    figure(2);
    imshow(Xhat2(:,:,11,10))
    figure(3);
    imshow(Y(:,:,11,10))
