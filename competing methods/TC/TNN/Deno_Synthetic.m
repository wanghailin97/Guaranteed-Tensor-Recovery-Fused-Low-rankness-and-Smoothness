%%
clear all
addpath(genpath(cd));
randn('seed',2013);
randn('seed',2013);

n  = 50;
n1 = n;
n2 = n;
n3 = n;
n4 = n;
r = floor(n*0.2); % t-SVD rank



U{1}=sqrt(n3)*dctmtx(n3);
U{2}=sqrt(n4)*dctmtx(n4);
A=randn(n1,r,n3,n4)/n;
B=randn(r,n2,n3,n4)/n;
XT =  htprod_U(A,B,U);
%XT =  htprod_fft(A,B);

[n1, n2, n3, n4]=size(XT);
Nways=size(XT);



sr = 0.72;
fprintf('Sampling ratio = %0.8e\n',sr);
temp = randperm(prod(Nways));
kks = round((sr)*prod(Nways));
mark = zeros((Nways)); 
mark(temp(1:kks)) = 1;
ZZ_=XT.*mark;



opts.mu = 1e-3;
opts.max_mu = 1e10;
opts.tol = 1e-8; 
opts.rho = 1.1;
opts.DEBUG = 1;
opts.max_iter =500;

%%  
  fprintf('===== t-SVD by Discrete Cosine Transform =====\n');
   tic
     [Xhat ,~,trank] = HTNN_U(U,XT,mark,opts);
   %[Xhat ,~,trank] = HTNN_FFT(XT,mark,opts);
   toc


RSE = norm(XT(:)-Xhat(:))/norm(XT(:));

fprintf('\nsampling rate: %f\n', sr);
fprintf('tubal rank of the underlying tensor: %d\n',r);
fprintf('tubal rank of the recovered tensor: %d\n', trank);
fprintf('relative recovery error: %.4e\n', RSE);
