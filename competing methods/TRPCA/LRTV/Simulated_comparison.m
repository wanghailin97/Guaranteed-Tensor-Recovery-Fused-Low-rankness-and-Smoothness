clear all;
% close all;
% clc;	


addpath ./PROPACK
addpath ./prox_operators 
%% simulated experiment
%----------------------------------------------load image---------------------------------------------------------
load Pavia_80.mat
load Pavia_noise.mat               % for case 2
% ratio = 0.2*ones(1,80);            % for case 1
% noiselevel = 0.1*ones(1,80);      % for case 1
% 
% load WDC.mat; OriData3 = WDC;
% load WDC_noise.mat                   % for case 2
% load stripes_WDC.mat
% ratio = 0.05*ones(1,191);            % for case 1
% noiselevel = 0.025*ones(1,191);      % for case 1
% %----------------------------------------------noise simulated---------------------------------------------------------
oriData3_noise = OriData3;
[M N p] = size(OriData3);
% Gaussian noise
for i =1:p
     oriData3_noise(:,:,i)=OriData3(:,:,i)  + noiselevel(i)*randn(M,N);
%        mixed(vv,:)=addnoise(mixed(vv,:),'mg',ratio);
end
% S&P noise
for i =1:p
     oriData3_noise(:,:,i)=imnoise(oriData3_noise(:,:,i),'salt & pepper',ratio(i));
end
% output_image = oriData3_noise;  
%% LRTV denoising
% tau = 0.005;%%0.005
% lambda =10/sqrt(M*N);
% rank =6;
% [ output_image out_value] = LRTV_accelerate(oriData3_noise, tau,lambda, rank);
%% local low rank combined global 3-D TV (LLRGTV)
par.lambda = 0.20;
par.tau  = 0.005;
par.r =2;
par.blocksize = 20;
par.stepsize  = 10;
% par.lambda = 5/par.blocksize;
% par.tau = 0.01;
% par.r = 4;
par.maxIter = 50;
par.tol = 1e-6;
[ output_image out_value] = LLRGTV(oriData3_noise,OriData3, par);

%% LRRSDS
% addpath code_LRRSDS
%         lambda = 0.1; % parameter for low rank regularization
%         lambda_s = 0.01; % parameter for sparse regularization
%         Desired_rank = 1;
%         AL_iter = 80;
% %figure;imshow(W,[]);
% [output_image, S, ~] = HSI_3D_LowRank_SpectralDifference(oriData3_noise,'MU',0.8, ...
%     'LAMBDA', lambda,'LAMBDA_S', lambda_s,'DESIRED_RANK',Desired_rank,'AL_ITERS',AL_iter, 'VERBOSE','yes','TRUE_X',oriData3_noise);
%% FastHyDe
% % estimate the location of sparse noise
% lambda =3.8/sqrt(M*N);
% tol = 1e-8;
% maxIter = 500;
% [~,sparsenoise] = inexact_alm_rpca(reshape(oriData3_noise,M*N,p), lambda, tol, maxIter);
% sparsenoise = reshape(sparsenoise,M,N,p);
% % oriData3_noise = reshape(oriData3_noise,M,N,p);
% idx = find(abs(sparsenoise)>=0.01);
% Mlocation = ones(M,N,p);
% Mlocation(idx) = 0;
% 
% noise_type='additive';
% iid = 0;
% addpath('./FastHyDe/FastHyDe');
% addpath('./FastHyDe/BM3D');
% p_subspace = 10;
% 
% [output_image, time_fasthyde] = FastHyIn(oriData3_noise,Mlocation, noise_type, iid, p_subspace);
% % [output_image, time_fasthyde] = FastHyDe(oriData3_noise,  noise_type, iid, p_subspace);
%% NAILRMR  
%   addpath ./NAIRSVD_public
%   addpath './NAIRSVD_public/demo_HySime'
% % parameter 
%  blocksize=20;
%  stepsize=8;
% %  opts.r = 3;
% %  opts.s = 0.15;
% %  opts.noiselevel = noiselevel';
% %  lambda = 1;
% % [oriData3_noise ] = reshape(inexact_alm_rpca(reshape(oriData3_noise,M*N,p)),M,N,p);
% % [output_image,~,~] =NAILRMR_denosing(oriData3_noise,OriData3,blocksize,stepsize,1,opts);  %NAIRSVD
% [output_image,~,~] = NAIRLMA_denosing(oriData3_noise,blocksize,stepsize,1);
% % 
% % [NAILRMA_PSNR,NAILRMA_SSIM,NAILRMA_SAM,NAILRMA_MQ] = evaluate(OriData3,output_image,M,N);
% % disp(['Method Name:NAILRMA   ', ', MPSNR=' num2str(mean(NAILRMA_PSNR),'%5.2f')  ...
% %       ',MSSIM = ' num2str(mean(NAILRMA_SSIM),'%5.4f')  ',SAM=' num2str(NAILRMA_SAM,'%5.2f')]);
%%
%    addpath 'NonLRMA for HSI denoising'
%     blocksize = 20; 
%     stepsize  = 4;
%     strname = 'Lap';
%     tic
%    [ output_image ] = NonLRMA_HSIdenoise( oriData3_noise,blocksize,stepsize,strname);
%     NonLRMA_time = toc;
%% 3-D video
% % [oriData3_noise ] = reshape(inexact_alm_rpca(reshape(oriData3_noise,M*N,p)),M,N,p);
% opts.beta    = [1 1 8];
% opts.print   = true;
% opts.method  = 'l2';
% 
% % Setup mu
% mu           =30;
% 
% % Main routine
% % tic
% output = deconvtvl2(oriData3_noise,1,mu,opts);output_image = output.f;
% % toc
%% ALM + SSAHTV
% lambda =1/sqrt(M*N);
% tol = 1e-8;
% maxIter = 500;
% [output_image] = inexact_alm_rpca(oriData3_noise, lambda, tol, maxIter); 
% 
% lambda = 2;
% sa=1;
% output_image = tvdenoise(output_image,lambda,sa,'l1');
%% SSAHTV  need tvrestore.m tvrestoresa.m shrink2.m uupdategs.m getWeights.m
% addpath './SSAHTV'
% % [oriData3_noise ] = reshape(inexact_alm_rpca(reshape(oriData3_noise,M*N,p)),M,N,p);
% lambda = 0.6;
% sa=1;
% output_image = tvdenoise(oriData3_noise,lambda,sa,'l1');
%% LRMR
% r = 4;
% slide =20;
% s= 0.10;
% % s = 0.08;
% stepsize = 4;
% tic
% [ output_image ] = LRMR_HSI_denoise( oriData3_noise,r,slide,s,stepsize );
% toc
%% tensor dictionary learning
% [oriData3_noise ] = reshape(inexact_alm_rpca(reshape(oriData3_noise,M*N,p)),M,N,p);
% addpath(genpath('./tensor_dl'));
% fprintf('Denoising by tensor dictionary learning ...\n');
% vstbmtf_params.peak_value = 1;
% vstbmtf_params.nsigma = mean(noiselevel);
% output_image = TensorDL(oriData3_noise, vstbmtf_params);

%% BM4D 
% addpath './BM4D_v2p2'
% % [oriData3_noise ] = reshape(inexact_alm_rpca(reshape(oriData3_noise,M*N,p)),M,N,p);
% tic
% [PSNR2, y_est, sigma_est] = bm4d(1, oriData3_noise, 0,0);
% toc 
% output_image = double(y_est);

%%
PSNRvector=zeros(1,p);
for i=1:1:p
    J=255*OriData3(:,:,i);

    I=255*output_image(:,:,i);

      PSNRvector(1,i)=PSNR(J,I,M,N);
end
% dlmwrite('PSNRvector.txt',PSNRvector,'delimiter','\t','newline','pc');
% 
PSNRvec = mean(PSNRvector)
SSIMvector=zeros(1,p);
for i=1:1:p
    J=255*OriData3(:,:,i);
%     Jnoise=oriData3_noise(:,:,i);
    I=255*output_image(:,:,i); 
%      [ SSIMnoisevector(1,i),ssim_map] = ssim(J,Jnoise);
      [ SSIMvector(1,i),ssim_map] = ssim(J,I);
end
mean(SSIMvector)
% dlmwrite('SSIMvector.txt',SSIMvector,'delimiter','\t','newline','pc');
