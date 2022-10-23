% This script run the TV based completion method

% You can:
%       1. Change test data by simply modifying variable 'tensor_num'
%       2. Change sampling rate by modifying variables  'sr '

% Please cite our paper listed in the following BibTex if you use any part of our source code or data:
%  (1)    
%      @article{JiTV,
%      title={Tensor completion using total variation and low-rank matrix factorization},
%      author={Ji, Teng-Yu and Huang, Ting-Zhu and Zhao, Xi-Le and Ma, Tian-Hui and Liu, Gang},
%      journal={Information Sciences},
%      volume={326},
%      pages={243--257},
%      year={2016},
%      publisher={Elsevier}
%      }
% 
%  (2)    
%      @article{Ji2017non,
%      title={A non-convex tensor rank approximation for tensor completion},
%      author={Ji, Teng-Yu and Huang, Ting-Zhu and Zhao, Xi-Le and Ma, Tian-Hui and Deng, Liang-Jian},
%      journal={Applied Mathematical Modelling},
%      volume={48},
%      pages={410--422},
%      year={2017},
%      publisher={Elsevier}
%      }
% 
%  (3)    
%      @article{JiangFramelet,
%      title={Matrix factorization for low-rank tensor completion using framelet prior},
%      author={Jiang, Tai-Xiang and Huang, Ting-Zhu and Zhao, Xi-Le and Ji, Teng-Yu and Deng, Liang-Jian},
%      journal={Information Sciences},
%      year={2017},
%      note={Accepted}
%      }
       
% created by Teng-Yu Ji 
% 1/13/2018
clc; clear; close all;

addpath(genpath('./lib'))
addpath('./testData')

for tensor_num=1
    switch tensor_num
        case 1 
            load video_suzie.mat
        case 2 
            load MRI.mat
        case 3 
            load HSI.mat
    end
%%
%% Set sampling ratio
sr = 0.2;

%% Initial data
X=X(:,:,1:10);

% normalized data
if max(X(:))>1
X=X/max(X(:));
end

Nway=[size(X,1), size(X,2), size(X,3)];
n1 = size(X,1); n2 = size(X,2); 
frames=size(X,3);
ratio = 0.005;
R=AdapN_Rank(X,ratio);
Y_tensorT = X;

p = round(sr*prod(Nway));
known = randsample(prod(Nway),p); data = Y_tensorT(known);
[known, id]= sort(known); data= data(id);
Y_tensor0= zeros(Nway); Y_tensor0(known)= data;
%imname=[num2str(tensor_num),'_tensor0'];

%% Initialization of the factor matrices X and A
for n = 1:3
    coNway(n) = prod(Nway)/Nway(n);
end
for i = 1:3
    Y0{i} = Unfold(Y_tensor0,Nway,i);
    Y0{i} = Y0{i}';
    X0{i}= rand(coNway(i), R(i));
    A0{i}= rand(R(i),Nway(i));
end
%    save(['./results/', imname],'Y_tensor0','Y_tensorT','known','A0','Y0','X0','Nway');
    
%% Initialization of the parameters
opts=[];
opts.maxit=2000;
opts.Ytr= Y_tensorT;
opts.tol=1e-4;
alpha=[1,1,1];
opts.alpha = alpha / sum(alpha);
rho=0.1;
opts.rho1=rho;
opts.rho2=rho;
opts.rho3=rho;
opts.mu=1; 
opts.beta=10;
opts.initer=10;
opts.miter=20;

%% Begin the comlpetion with MF_TV
        fprintf('\n');
        disp('Begin the comlpetion with TV')
        t0= tic;
        [Y_TV, A, X, Out]= LRTC_TV(Y0, data, A0, X0,Y_tensor0, Nway, known, opts, n1, n2);
        time= toc(t0);
        psnr_one=PSNR(Y_TV, Y_tensorT);
        imname=[num2str(tensor_num),'result_TV_psnr_',num2str(psnr_one), '.mat'];
%        save(['./results/', imname],'A','X','Y_TV','opts','Out');
end