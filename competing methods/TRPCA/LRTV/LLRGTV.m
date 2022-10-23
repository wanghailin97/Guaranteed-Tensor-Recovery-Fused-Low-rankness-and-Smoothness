function [ output_image,out_value] = LLRGTV(oriData3_noise,OriData3, par)
 % Solve problem
% solve the proposed LLRGTV problem
% -------------------------------------------------------------------------------------------------------------
% Reference paper: ¡°Hyperspectral Image Denoising Using Local Low-rank Matrix Recovery and Global Total Variation¡± 
% Author: Wei He (Oct., 2016)
% E-mail addresses:(weihe1990@whu.edu.cn)
% --------------------------------------------------INPUT-----------------------------------------------------
% par.lambda                                sparse regularization parameter (recommend 0.2) 
% par.tau                                   SSTV regularization parameter(recommend 0.005)
% par.r                                     upper bound rank (estimated by multiple regression theory-based approach)
% par.blocksize                             size of each small cube (20*20)
% par.stepsize  = 8;                        stepsize
% par.maxIter = 20;                         max iteration
% par.tol = 1e-6;                           convergence conditions
% --------------------------------------------------OUTPUT-----------------------------------------------------
%  output_iamge                              3-D denoised image
%  out_value                                 MPSNR and MSSIM valuses of each iteration 
% -------------------------------------------------------------------------------------------------------------
% Lp Jp and Sp are the tensors collecting the patches of 3D L, J and
% S,respectively
% -------------------------------------------------------------------------------------------------------------
blocksize = par.blocksize;
stepsize  = par.stepsize;
lambda    = par.lambda;
tau       = par.tau;
r         = par.r;
maxIter   = par.maxIter;
tol       = par.tol;
% the parameters in this part is the link between global and local iamge
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 
[M,N,p] = size(oriData3_noise);

R         =   M-blocksize+1;
C         =   N-blocksize+1;
rr        =   [1:stepsize:R];
% rr        =   [rr rr(end)+1:R];
if rr(end)<R;rr = [rr R];else end;
cc        =   [1:stepsize:C];
% cc        =   [cc cc(end)+1:C];
if cc(end)<R;cc = [cc C];else end;
row       =   length(rr);
column    =   length(cc);
% Index image
Idxx     =   (1:row*column);
Idxx     =   reshape(Idxx, row, column);

par.rr   = rr;
par.cc   = cc;
par.row  = row;
par.column = column;
par.Idxx   = Idxx;
par.imagesize = [M,N,p];
%% Initializing optimization variables
Op = imageTopatch(oriData3_noise,par); % store 3-D image patch by patch
weight = calweight(par);

[M1,N1,p1] = size(Op);% M1 = stepsize.^2, N1 = p; p1 number of patches; 

Lp = randn(M1,N1,p1); Jp = Lp; Sp = zeros(M1,N1,p1); X = patchToimage(Jp,par);
YO = zeros(M1,N1,p1); YL = zeros(M1,N1,p1); YX = zeros(M,N,p);

u1 = zeros(M,N,p); u2 = zeros(M,N,p); u3 = zeros(M,N,p);
y1 = zeros(M,N,p); y2 = zeros(M,N,p); y3 = zeros(M,N,p);

% for the 3D-TV norm
% define operators
% in this problem H=1;beta = [1 1 0.5];
H =1;beta = [1 1 0.5];  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
%% 
rho = 1.5;
max_mu1 = 1e6;
mu = 1e-2;
sv =10;
%% 
out_value = [];
out_value.SSIM =[];
out_value.PSNR =[];
out_value.obj= [];
%% main loop
iter = 0;
tic
while iter<maxIter
    iter = iter + 1;
    
%     Update L & Lp
temp = (Op - Sp + Jp)/2 + (YO - YL) /(2*mu);
for i =1:p1
%     Update Lp
  temp1 = temp(:,:,i);
  if  choosvd(p,sv) ==1
      [U, sigma, V] = lansvd(temp1, sv, 'L');
  else
      [U,sigma,V] = svd(temp1,'econ');
  end
    sigma = diag(sigma);
    svp = min(length(find(sigma>1/(2*mu))),r);
    if svp<sv
        sv = min(svp + 1, p);
    else
       sv = min(svp + round(0.05*p), p);  
    end
  Lp(:,:,i) = U(:, 1:svp) * diag(sigma(1:svp) - 1/(2*mu)) * V(:, 1:svp)'; 

end

%    Update S & Sp
  temp_S = Op - Lp+ YO/mu;
  Sp_hat = max(temp_S - lambda/mu, 0);
  Sp = Sp_hat+min(temp_S + lambda/mu, 0);

%    Update J & Jp
 temp_C = Lp + YL/mu;
 temp_D = X - YX/mu;
 J      = update_M(temp_C,temp_D,weight,par);
 Jp     = imageTopatch(J,par);

%       Update X
temp = J + YX/mu;
Htg         = imfilter(temp, H, 'circular');
rhs   = fftn(Htg + Dt(u1+(1/mu)*y1,  u2+(1/mu)*y2, u3+(1/mu)*y3));
eigA  = (mu/mu)*eigHtH + eigDtD + eigEtE;
X     = real(ifftn(rhs./eigA));
% X = reshape(f,M*N,p);

% update U = [u1 u2 u3]
 [Df1 Df2 Df3] = D(X);

  v1 = Df1-(1/mu)*y1;
  v2 = Df2-(1/mu)*y2;
  v3 = Df3-(1/mu)*y3;
  u1 = max(v1 - tau/mu, 0); u1 = u1 + min(v1 + tau/mu,0);
  u2 = max(v2 - tau/mu, 0); u2 = u2 + min(v2 + tau/mu,0);
  u3 = max(v3 - tau/mu, 0); u3 = u3 + min(v3 + tau/mu,0);

  leq1 = Op - Lp - Sp;
  leq2 = Lp - Jp;
  leq3 = J - X;
%% stop criterion          
%     stopC = max(max(max(abs(leq1))),max(max(abs(leq2))));
 out_value.obj1(iter) = sum(abs(leq1(:)));
 out_value.obj2(iter) = sum(abs(leq2(:)));
 out_value.obj3(iter) = sum(abs(leq3(:)));
stopC1 = max(abs(leq1(:)));
stopC2 = max(abs(leq2(:)));
stopC3 = max(abs(leq3(:)));
stopC = max(max(stopC1,stopC2),stopC3);

% disp(['iter ' num2str(iter) ',mu=' num2str(mu,'%2.1e')  ...
%            ',stopALM=' num2str(stopC,'%2.3e')]);
 
 if stopC<tol
   break;
 else
   YO = YO + mu*leq1;
   YL = YL + mu*leq2;
   YX = YX + mu*leq3;
 
   y1   = y1 + mu*(u1 - Df1);
   y2   = y2 + mu*(u2 - Df2);
   y3   = y3 + mu*(u3 - Df3);
   mu = min(max_mu1,mu*rho);
   
 end
%% output SSIM and PSNR values of each step

output_image = J;
   PSNRvector=zeros(1,p);
for i=1:1:p
    J=255*OriData3(:,:,i);

    I=255*output_image(:,:,i);

      PSNRvector(1,i)=PSNR(J,I,M,N);
end
% dlmwrite('PSNRvector.txt',PSNRvector,'delimiter','\t','newline','pc');
% 
out_value.PSNR(iter) = mean(PSNRvector);
SSIMvector=zeros(1,p);
for i=1:1:p
    J=255*OriData3(:,:,i);
%     Jnoise=oriData3_noise(:,:,i);
    I=255*output_image(:,:,i); 
%      [ SSIMnoisevector(1,i),ssim_map] = ssim(J,Jnoise);
      [ SSIMvector(1,i),ssim_map] = ssim(J,I);
end
out_value.SSIM(iter)=mean(SSIMvector); 
%% output SSIM and PSNR values of each step  
% output_image = J;
% 
% PSNRvector=zeros(1,p);
% for i=1:1:p
%     J=255*OriData3(:,:,i);
% 
%     I=255*output_image(:,:,i);
% 
%       PSNRvector(1,i)=PSNR(J,I,M,N);
% end
% % dlmwrite('PSNRvector.txt',PSNRvector,'delimiter','\t','newline','pc');
% % 
% PSNRvec = mean(PSNRvector)
% SSIMvector=zeros(1,p);
% for i=1:1:p
%     J=255*OriData3(:,:,i);
% %     Jnoise=oriData3_noise(:,:,i);
%     I=255*output_image(:,:,i); 
% %      [ SSIMnoisevector(1,i),ssim_map] = ssim(J,Jnoise);
%       [ SSIMvector(1,i),ssim_map] = ssim(J,I);
% end
% mean(SSIMvector)
end
toc

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

function [ weight] = calweight(par )
% this function is to calculate the weight from 3-D image to blockes
% by Wei He 
%%
% [M1 N1 p1] = size(outpatch);
blocksize  = par.blocksize;
imagesize  = par.imagesize;

rr        = par.rr;
cc        = par.cc;
row       = par.row;
column    = par.column;
Idxx      = par.Idxx;
Numofpatch= row * column; 

weight   = zeros(imagesize);
for idx = 1:Numofpatch
    [rowidx,columnidx] = find(Idxx==idx); % find the location in index image
    i = rr(rowidx); j = cc(columnidx);    % find the location in original image
    weight(i:i+blocksize-1,j:j+blocksize-1,:)= weight(i:i+blocksize-1,j:j+blocksize-1,:)+ 1;
end  
end