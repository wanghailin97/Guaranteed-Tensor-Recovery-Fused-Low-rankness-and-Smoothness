function [ output_image,out_value] = LRTV_denoising(oriData3_noise, tau,lambda,G1,G0, r)
 % Solve problem
% solve the following problem (TV regularized LR and MC problem)
%  argmin   ||L||_nuclear + tao * ||X||_TV+lanbda*||E||_1
%                       s.t. X = L L = TM and P_omega(D) = Pomega_(TM + E) and rank(L)<=r;
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
D = zeros(M*N,p);
for i=1:p 
  bandp  = oriData3_noise(:,:,i);
  D(:,i) = bandp(:); 
end
[d p] = size(D);
%%%%%% 计算D在非零位置的F范数
d_norm = norm(D, 'fro');
%%%%%% 计算D在非零位置的F范数
% initialize
% b_norm = norm(D, 'fro');
tol = 1e-6;
tol1 = tol;
tol2 = tol1;
maxIter = 40;

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
X = L;TM = X;
E = sparse(d,p);

Y1 = zeros(d,p);
Y2 = zeros(d,p);
Y3 = zeros(d,p);


% for the TV norm
param2.verbose=1;
param2.max_iter=10;
param2.verbose = 0;
% g1.prox=@(x, T) prox_TV(x, T*mu1/(2*tau), param2);
% g1.norm=@(x) tau*TV_norm(x); 


% main loop
iter = 0;
tic
while iter<maxIter
    iter = iter + 1;      
%     Updata L
temp = (mu1*X +mu2* TM - (Y1+Y2))/(mu1+mu2);
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

%       Updata X
temp = L + Y1/mu1;
X = L;
for i =1:p  
 z = prox_TV(reshape(temp(:,i),[M,N]),2*tau/mu1,param2);
%    z=denoise_bound(reshape(temp(:,i),[M,N]),tau/mu1,0,1);
%   z = prox_TV(reshape(temp(:,i),[M,N]),tau/mu1,param2);
 X(:,i) = z(:);  
end

%   Update TM
mid2 = (L+Y2/mu2);
mid3 = D-E+Y3/mu3;
TM(G0) = mid2(G0);
TM(G1) = (mid2(G1) + mid3(G1))/2;

% updata E
temp_E = D - TM + Y3/mu3;
temp_E(G0) = 0;
E_hat = max(temp_E - lambda/mu3, 0);
E = E_hat+min(temp_E + lambda/mu3, 0);

leq1 = L-X;
leq2 = L-TM;
leq3 = D -TM -E ;
leq3(G0) = 0;
%% stop criterion          
%     stopC = max(max(max(abs(leq1))),max(max(abs(leq2))));
 
stopC1 = max(max(abs(leq1)));
stopC2 = max(max(abs(leq1)));
stopC3 = norm(leq3, 'fro') / d_norm;
disp(['iter ' num2str(iter) ',mu=1' num2str(mu1,'%2.1e')  ...
           ',rank = ' num2str(rank(L))  ',stopALM=' num2str(stopC3,'%2.3e')...
           ',stopE=' num2str(max(stopC1,stopC2),'%2.3e')]);
 
 if stopC2<tol  & stopC1<tol1 & stopC3<tol2
   break;
 else
   Y1 = Y1 + mu1*leq1;
   Y2 = Y2 + mu2*leq2;
   Y3 = Y3 + mu3*leq3;
   mu1 = min(max_mu1,mu1*rho);
   mu2 = min(max_mu1,mu2*rho);
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
end

