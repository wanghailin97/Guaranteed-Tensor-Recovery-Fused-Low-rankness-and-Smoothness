function [ output_image output_sparse output_noise] = LRTV_accelerate(oriData3_noise, tau,lambda, r)
 % Solve problem
% solve the following problem (TV regularized LR and MC problem)
%  argmin   ||L||_nuclear + tao * ||X||_TV+lanbda*||E||_1
%                       s.t. X = L D = L + E and rank(L)<=r;
%     via IALM 
% -------------------------------------------------------------------------------------------------------------
% Reference paper: W. He, H. Zhang, L. Zhang, and  H. Shen, “Total-Variation-Regularized
% Low-Rank Matrix Factorization for Hyperspectral Image Restoration,” IEEE Trans. Geosci. Remote Sens., 
% vol. 54, pp. 178-188, Jan. 2016.
% Author: Wei He (November, 2014)
% E-mail addresses:(weihe1990@whu.edu.cn)
% --------------------------------------------------INPUT-----------------------------------------------------
%  oriData3_noise                            noisy 3-D image of size M*N*p normalized to [0,1] band by band
%  tao                                       (recommended value 0.01)              
%  lambda 
%  G1(omega)                                 all one matrix of size (M*N)*p
%  G0(omega~)                                all zero matrix of size (M*N)*p
%  r                                         rank constraint
% --------------------------------------------------OUTPUT-----------------------------------------------------
%  output_iamge                              3-D denoised image
%  out_value                                 MPSNR and MSSIM valuses of each iteration 
% -------------------------------------------------------------------------------------------------------------
% Note: the parameters G0 and G1 do not have any effect. these two
% parameters are used to solve the impainting problem with the location of
% missing pixels to be known.
% -------------------------------------------------------------------------------------------------------------
[M N p] = size(oriData3_noise);
Dob = zeros(M*N,p);
for i=1:p 
  bandp  = oriData3_noise(:,:,i);
  Dob(:,i) = bandp(:); 
end
[d p] = size(Dob);
%%%%%% 计算D在非零位置的F范数
d_norm = norm(Dob, 'fro');
%%%%%% 计算D在非零位置的F范数
% initialize
% b_norm = norm(D, 'fro');
tol = 1e-6;
tol1 = tol;
tol2 = tol1;
maxIter = 30;

rho = 1.5;
max_mu1 = 1e6;
mu1 = 1e-2;
mu2  = mu1;
mu3 = mu2;
sv =10;
%% 
out_value = [];
out_value.SSIM =[];
out_value.PSNR =[];
%% Initializing optimization variables
% intialize
% L = zeros(d,p);
% X = zeros(d,p);
% TM = zeros(d,p);
L = rand(d,p);
X = L;
E = sparse(d,p);

Y1 = zeros(d,p);
Y2 = zeros(d,p);

u1 = zeros(M,N,p); u2 = zeros(M,N,p); u3 = zeros(M,N,p);
y1 = zeros(M,N,p); y2 = zeros(M,N,p); y3 = zeros(M,N,p);

% for the 3D-TV norm
% define operators
% in this problem H=1;beta = [1 1 10];
H =1;beta = [1 1 0];
eigHtH      = abs(fftn(H, [M N p])).^2;
eigDtD      = abs(beta(1)*fftn([1 -1],  [M N p])).^2 + abs(beta(2)*fftn([1 -1]', [M N p])).^2;
if p>1
    d_tmp(1,1,1)= 1; d_tmp(1,1,2)= -1;
    eigEtE  = abs(beta(3)*fftn(d_tmp, [M N p])).^2;
else
    eigEtE = 0;
end
% Htg         = imfilter(g, H, 'circular');
[D,Dt]      = defDDt(beta);

% main loop
iter = 0;
tic
while iter<maxIter
    iter = iter + 1;      
%     Update L
temp = (mu1*X +mu2* (Dob-E) + (Y1+Y2))/(mu1+mu2);
if  choosvd(p,sv) ==1
    [U, sigma, V] = lansvd(temp, sv, 'L');
else
    [U,sigma,V] = svd(temp,'econ');
end
    sigma = diag(sigma);
    svp = min(length(find(sigma>1/(mu1+mu2))),r);
    if svp<sv
        sv = min(svp + 1, p);
    else
       sv = min(svp + round(0.05*p), p);
        
    end
 L = U(:, 1:svp) * diag(sigma(1:svp) - 1/(mu1+mu2)) * V(:, 1:svp)'; 

 
%       Update X
temp = reshape(L - Y1/mu1,M,N,p);
Htg         = imfilter(temp, H, 'circular');
rhs   = fftn((mu1/mu3)*Htg + Dt(u1+(1/mu3)*y1,  u2+(1/mu3)*y2, u3+(1/mu3)*y3));
eigA  = (mu1/mu3)*eigHtH + eigDtD + eigEtE;
f     = real(ifftn(rhs./eigA));
X = reshape(f,M*N,p);

% update M = [u1 u2 u3]
 [Df1 Df2 Df3] = D(f);

%   v1 = Df1-(1/mu3)*y1;
%   v2 = Df2-(1/mu3)*y2;
%   v3 = Df3-(1/mu3)*y3;
%   v  = sqrt(v1.^2 + v2.^2 + v3.^2);
%   v(v==0) = 1e-6;
%   v  = max(v - tau/mu3, 0)./v;
%   u1 = v1.*v;
%   u2 = v2.*v;
%   u3 = v3.*v;

  v1 = Df1-(1/mu3)*y1;
  v2 = Df2-(1/mu3)*y2;
  v3 = Df3-(1/mu3)*y3;
  u1 = max(v1 - tau/mu3, 0); u1 = u1 + min(v1 + tau/mu3,0);
  u2 = max(v2 - tau/mu3, 0); u2 = u2 + min(v2 + tau/mu3,0);
  u3 = max(v3 - tau/mu3, 0); u3 = u3 + min(v3 + tau/mu3,0);

% update E
temp_E = Dob - L + Y2/mu2;
E_hat = max(temp_E - lambda/mu2, 0);
E = E_hat+min(temp_E + lambda/mu2, 0);

leq1 = X - L;
leq2 = Dob -L -E ;
%% stop criterion          
%     stopC = max(max(max(abs(leq1))),max(max(abs(leq2))));
 
stopC1 = max(max(abs(leq1)));
stopC2 = norm(leq2, 'fro') / d_norm;
disp(['iter ' num2str(iter) ',mu=' num2str(mu1,'%2.1e')  ...
           ',rank = ' num2str(rank(L))  ',stopALM=' num2str(stopC2,'%2.3e')...
           ',stopE=' num2str(stopC1,'%2.3e')]);
 
 if stopC1<tol  & stopC2<tol2
   break;
 else
   Y1 = Y1 + mu1*leq1;
   Y2 = Y2 + mu2*leq2;
   mu1 = min(max_mu1,mu1*rho);
   mu2 = min(max_mu1,mu2*rho); 
   
   y1   = y1 + mu3*(u1 - Df1);
   y2   = y2 + mu3*(u2 - Df2);
   y3   = y3 + mu3*(u3 - Df3);
   mu3 = min(max_mu1,mu3*rho);
   
 end
%% output SSIM and PSNR values of each step
% load simu_indian
% OriData3 = simu_indian;
% 
% output_image = reshape(L,[M,N,p]);
%    PSNRvector=zeros(1,p);
% for i=1:1:p
%     J=255*OriData3(:,:,i);
% 
%     I=255*output_image(:,:,i);
% 
%       PSNRvector(1,i)=PSNR(J,I,M,N);
% end
% % dlmwrite('PSNRvector.txt',PSNRvector,'delimiter','\t','newline','pc');
% % 
% out_value.PSNR(iter) = mean(PSNRvector);
% SSIMvector=zeros(1,p);
% for i=1:1:p
%     J=255*OriData3(:,:,i);
% %     Jnoise=oriData3_noise(:,:,i);
%     I=255*output_image(:,:,i); 
% %      [ SSIMnoisevector(1,i),ssim_map] = ssim(J,Jnoise);
%       [ SSIMvector(1,i),ssim_map] = ssim(J,I);
% end
% out_value.SSIM(iter)=mean(SSIMvector); 
%% output SSIM and PSNR values of each step    
end
toc
output_image = reshape(L,[M,N,p]);
output_sparse = reshape(E,[M,N,p]);
output_noise = oriData3_noise-output_image-output_sparse;
end

function [D,Dt] = defDDt(beta)
D  = @(U) ForwardD(U, beta);
Dt = @(X,Y,Z) Dive(X,Y,Z, beta);
end

function [Dux,Duy,Duz] = ForwardD(U, beta)
frames = size(U, 3);
Dux = beta(1)*[diff(U,1,2), U(:,1,:) - U(:,end,:)];
Duy = beta(2)*[diff(U,1,1); U(1,:,:) - U(end,:,:)];
Duz(:,:,1:frames-1) = beta(3)*diff(U,1,3); 
Duz(:,:,frames)     = beta(3)*(U(:,:,1) - U(:,:,end));
end

function DtXYZ = Dive(X,Y,Z, beta)
frames = size(X, 3);
DtXYZ = [X(:,end,:) - X(:, 1,:), -diff(X,1,2)];
DtXYZ = beta(1)*DtXYZ + beta(2)*[Y(end,:,:) - Y(1, :,:); -diff(Y,1,1)];
Tmp(:,:,1) = Z(:,:,end) - Z(:,:,1);
Tmp(:,:,2:frames) = -diff(Z,1,3);
DtXYZ = DtXYZ + beta(3)*Tmp;
end

