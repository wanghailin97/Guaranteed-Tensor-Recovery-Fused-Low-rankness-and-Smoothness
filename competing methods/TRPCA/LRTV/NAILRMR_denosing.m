function [output_image,rsvdPSNR,rsvdSSIM] = NAILRMR_denosing(oriData3_noise,OriData3,blocksize,stepsize,power,opts)
%   reference paper He. Wei, Zhang. Hongyan, Zhang. Liangpei, and Shen. Huanfeng,
%   "Hyperspectral Image Denoising via Noise-Adjusted Iterative Low-Rank Matrix Approximation," 
%   J-star, vol. 8, pp. 3050-3061, 2015.
% addpath 'demo_HySime'
%
% oriData3_noise          输入噪声影像
% OriData3                参考影像（仅用于计算PSNR值以及SSIM值）
% blocksize               子块的大小
% stepsize                步长
% r                       每个子块的维数估计
% delta                   迭代因子估计
% power                   set to 1

rsvdPSNR = [];
rsvdSSIM = [];
M = blocksize; % the size of each block
stepszie=stepsize;
[m,n,p] = size(oriData3_noise); % the size of image
noiseimage = oriData3_noise;
% 估计噪声 输出秩的最大化估计以及迭代因子
[m,n,p] = size(oriData3_noise);
esty = zeros(m*n,p);
for i= 1:p
    mid = oriData3_noise(:,:,i);
    esty(:,i) = mid(:);
end
esty =esty';
noise_type = 'additive';
verbose ='off';
[w Rn] = estNoise(esty,noise_type,verbose);
r = opts.r;
s = opts.s;
delta = exp(-5.*opts.noiselevel);

% rank estimation
% [Uesty,Sesty,Vesty] = svd(esty,'econ');
% [Uw,Sw,Vw]          = svd(w,'econ');
% diagSesty = diag(Sesty);
% r = length(find(diagSesty>=Sw(1,1))); 
% % opts.noiselevel
% r =7 ;s = 0.15;
% % noise adjusted parameter estimation
% Rn = sqrt(Rn);
% imageband = size(Rn,2);
% delta = zeros(1,imageband);
% diagRn = diag(Rn)';
% mu =5; delta = exp(-mu.*diagRn);
% % r =3 ;s = 0.1; diagRn = 0.1*ones(1,p);
% % mu =5; delta = exp(-mu.*diagRn);
% % 估计噪声 输出秩的最大化估计以及迭代因子       
converged = 0 ;
iter = 1;

%% slide the 3-D image into block by block 
R         =   m-M+1;
C         =   n-M+1;
rr        =   [1:stepszie:R];
rr        =   [rr rr(end)+1:R];
cc        =   [1:stepszie:C];
cc        =   [cc cc(end)+1:C];
row       =   length(rr);
column    =   length(cc);

tic
while ~converged
       
    if iter ==1 % the first iteration 
         oriData3_noise = oriData3_noise;        
    else
     for i =1:p % noise adjusted iteration
       oriData3_noise(:,:,i)  =   output_image(:,:,i) + delta(1,i)*(oriData3_noise(:,:,i) - output_image(:,:,i));
     end  
    end

clear_image=zeros(m,n,p);
Weight = zeros(m,n,p);
patch_block = zeros(M^2,p);

%%
  for   rownumber =1:row
     for columnnumber = 1:column
         i = rr(rownumber);
         j = cc(columnnumber);
        for  k=1:1:p   % reshape the patch to 2-D matrix,and calculate the relative weight matrix                
         patch_reference = oriData3_noise(i:i+M-1,j:j+M-1,k); 
         patch_block(:,k) =  patch_reference(:);
         Weight(i:1:i+M-1,j:1:j+M-1,k) = Weight(i:1:i+M-1,j:1:j+M-1,k)+1;             
        end
%         [clear_patch_block] = BRP(patch_block,r,power);% RSVD denoising patch by patch 
          [clear_patch_block,S,~,~]=SSGoDec( patch_block,r,s,power);
        for m2=1:1:p % put denoised patch back to the image
           clear_image(i:1:i+M-1,j:1:j+M-1,m2) = reshape(clear_patch_block(:,m2),M,M)+clear_image(i:1:i+M-1,j:1:j+M-1,m2);      
        end 
    end
  end
 Weight_last = 1./Weight;
 output_image= Weight_last.*clear_image;  % the last output image

%% Assessment
 PSNRvector=zeros(1,p);
for i=1:1:p
    J=255*OriData3(:,:,i);

    I=255*output_image(:,:,i);

      PSNRvector(1,i)=PSNR(J,I,m,n);
end

PSNRvec = mean(PSNRvector);
SSIMvector=zeros(1,p);
for i=1:1:p
    J=255*OriData3(:,:,i);
%     Jnoise=oriData3_noise(:,:,i);
    I=255*output_image(:,:,i); 
%      [ SSIMnoisevector(1,i),ssim_map] = ssim(J,Jnoise);
      [ SSIMvector(1,i),ssim_map] = ssim(J,I);
end
SSIMvec = mean(SSIMvector);

rsvdPSNR(1,iter) = PSNRvec;
rsvdSSIM(1,iter) = SSIMvec;

rel = norm(output_image(:) - oriData3_noise(:),2)/norm(oriData3_noise(:),2);
 if rel<=1e-3
   converged = true;
 end  
 if iter >=5
  converged = true;
 end
fprintf( 'iteration = %d,PSNR = %2.2f,SSIM = %2.4f,rel = %2.4f\n',iter,  PSNRvec,SSIMvec,rel );
iter = iter+1;
end
toc
end  

%%
function[L,S,RMSE,error]=SSGoDec(X,rank,card,power)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Semi-Soft GoDec Algotithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
%X: nxp data matrix with n samples and p features
%rank: rank(L)<=rank
%card: card(S)<=card
%power: >=0, power scheme modification, increasing it lead to better
%accuracy and more time cost
%OUTPUTS:
%L:Low-rank part
%S:Sparse part
%RMSE: error
%error: ||X-L-S||/||X||
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REFERENCE:
%Tianyi Zhou and Dacheng Tao, "GoDec: Randomized Lo-rank & Sparse Matrix
%Decomposition in Noisy Case", ICML 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tianyi Zhou, 2011, All rights reserved.

%iteration parameters
iter_max=1e+2;
error_bound=1e-10;
iter=1;
RMSE=[];

%matrix size
[m,n]=size(X);
if m<n X=X'; end

%initialization of L and S
L=X;
S=sparse(zeros(size(X)));

% tic;
while true

    %Update of L
    Y2=randn(n,rank);
    for i=1:power+1
        Y1=L*Y2;
        Y2=L'*Y1;
    end
    [Q,R]=qr(Y2,0);
    L_new=(L*Q)*Q';
    
    %Update of S
    T=L-L_new+S;
    L=L_new;
    S=wthresh(T,'s',card);
%     [~,idx]=sort(abs(T(:)),'descend');
%     S=zeros(size(X));S(idx(1:card))=T(idx(1:card));
    
    %Error, stopping criteria
    %T(idx(1:card))=0;
    T=T-S;
    RMSE=[RMSE norm(T(:))];
    if RMSE(end)<error_bound || iter>iter_max
        break;
    else
        L=L+T;
    end
    iter=iter+1;
    
end
% toc;

LS=L+S;
error=norm(LS(:)-X(:))/norm(X(:));
if m<n 
    LS=LS';
    L=L';
    S=S';
end
end