
clear; 
close all; 
% 
% addpath('prox_operators');
% addpath('dataset');
% addpath('quality_assess');
% addpath('tensor_toolbox');
% addpath('my_main');

%% data claire
% load('Cl')
% Ltr = Cl./max(Cl(:));
% Ltr =Ltr(:,:,1:10);
video=cell2mat(struct2cell(load('D:\目前工作\原数据\video6.mat')));

page = 25;

Ltr = video(:,:,:,1:page);
n = size(Ltr);
Ltr = reshape(Ltr,[n(1),n(2),n(3)*n(4)]);
%%
[n1,n2,n3] = size(Ltr);
 sr = 0.6;
 level = 0.5;
 sigma = 0.01;

Bn = imnoise(Ltr,'salt & pepper',level);
noisy = Bn + sigma*randn([n1,n2,n3]);

%% sampling ratio and Omega


fprintf('Sampling ratio = %0.8e\n',sr);
temp = randperm(n1*n2*n3);
kks = round((sr)*n1*n2*n3);
Omega = zeros(n1,n2,n3); 
Omega(temp(1:kks)) = 1;

Y = noisy.*Omega;

beta = 1;
gamma = 0.001;
lambda = 3.2/sqrt(sr*n3*max(n1,n2));
rho1 = 400;
tau = 1.618;
tol = 5e-3;
itemax = 50;
theta = 50;
xi = 1.2;

opts.theta  = theta;
opts.xi     = xi;
opts.tau    = tau;
opts.rho    = rho1;
opts.lambda = lambda;
opts.beta   = beta;
opts.gamma  = gamma;
opts.tol    = tol;   
opts.itemax = itemax;


%% main loop FFT
% fprintf('FFT \n');
% tic;
% [Lf,S,k,eta,beta] = LowSparN(Y, Omega, opts);
% toc;
%  %% Printf PSNR and SSIM
%     Error = norm(Lf(:)-Ltr(:))/norm(Lf(:));
%     fprintf('Relative error = %0.8e\n',Error);
%     PSNR = psnr(Ltr(:),Lf(:));
%     fprintf('PSNR = %0.8e\n',PSNR);
%     SSIM = ssim3d(Lf*255,Ltr*255);
%     fprintf('SSIM = %0.8e\n',SSIM);
%     
%     t=2;
% figure(1)
% subplot(2,2,1);imshow(Ltr(:,:,t),'Border','tight');
% subplot(2,2,2);imshow(Y(:,:,t),'Border','tight');
% subplot(2,2,3);imshow(S(:,:,t),'Border','tight');
% subplot(2,2,4);imshow(Lf(:,:,t),'Border','tight');
% DCT
 fprintf('DCT \n');

%% paraters            
beta = 1;
gamma = 0.1;
lambda = 80/sqrt(sr*n3*max(n1,n2));
rho2 = 1000;
opts.beta   = beta;
opts.gamma  = gamma;
opts.lambda = lambda;
opts.rho    = rho2;
% 
tic;
[Ld Sd kd etad betad] = LowSparNDCT(Y, Omega, opts);
toc;

% Printf PSNR and SSIM
    Error = norm(Ld(:)-Ltr(:))/norm(Ld(:));
    fprintf('Relative error = %0.8e\n',Error);
    PSNR = psnr(Ltr(:),Ld(:));
    fprintf('PSNR = %0.8e\n',PSNR);
    SSIM = ssim3d(Ld*255,Ltr*255);
    fprintf('SSIM = %0.8e\n',SSIM);
    
    
fprintf('Unitary transform\n');
%% U matrix

 rho = 1000;
 opts.rho    = rho;
O = Unfold(Ld,[n1 n2 n3],3);
[U D V] = svd(O,'econ');

clear D  V

tic;
[Lu Su ku etau betau] = LowSparNU(U, Y, Omega, opts);
toc;
etau
% Printf PSNR and SSIM
    Error = norm(Lu(:)-Ltr(:))/norm(Lu(:));
    fprintf('Relative error = %0.8e\n',Error);
    PSNR = psnr(Ltr(:),Lu(:));
    fprintf('PSNR = %0.8e\n',PSNR);
    SSIM = ssim3d(Lu*255,Ltr*255);
    fprintf('SSIM = %0.8e\n',SSIM);


