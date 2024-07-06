%% Demo of Tensor Completion
%  test for TCTV based RGB image inpainting and comparison with other SOTA methods
%==========================================================================
% This script compares many tensor completion methods
% listed as follows:
%   1. SNN
%   2. BCPF
%   3. KBR
%   4. IRTNN
%   5. TNN
%   6. MF+TV
%   7. SNN+TV
%   8. SPC+TV
%   9. TNN+TV
%   10. TCTV (proposed method)
%
% -- three quality assessment (QA) indices -- PSNR, SSIM, FSIM, 
% -- are calculated for each method's denoising result, see Img_QA
%
% You can:
%       1. Type 'Demo' to to run various methods and see the pre-computed results.
%       2. Change test Img by simply modifying variable 'dataName' in Demo.m 
%          (NOTE: make sure your Img meets the format requirements).
%       3. Change missing level by modifying variables  'missing_rate' in Demo.m
%       4. Select competing methods by turn on/off the enable-bits in Demo.m
%
% More detail can be found in 
% Hailin Wang et.al, Guaranteed Tensor Recovery Fused Low-rankness and Smoothness, 2022
%
% by Hailin Wang, 2022
%==========================================================================
% The inpainting results of the RGB image "Lena" when missing_rate = 90%
% ================== QA Results =====================
%    Method    MPSNR    MSSIM    MFSIM      Time  
%  Observed   5.578     0.027    0.498    0.000   
%       SNN   19.813    0.475    0.749    3.217   
%      BCPF   20.356    0.404    0.717    66.595   
%       KBR   18.756    0.327    0.674    28.062   
%     IRTNN   18.855    0.342    0.690    54.228   
%       TNN   19.721    0.436    0.734    2.970   
%     MF_TV   9.130     0.115    0.475    290.494   
%    SNN_TV   22.759    0.703    0.812    15.478   
%    SPC_TV   20.260    0.481    0.755    29.072   
%    TNN_TV   23.747    0.737    0.822    13.627   
%      TCTV   26.227    0.803    0.892    18.305 
% =========================================================================
% The inpainting results of the RGB image "Einstein" when missing_rate = 99%
% ================== QA Results =====================
%    Method    MPSNR    MSSIM    MFSIM      Time  
%  Observed   11.294    0.038    0.598    0.000   
%       SNN   15.852    0.408    0.709    30.437   
%      BCPF      -        -        -         -    % memory overflow
%       KBR   16.459    0.338    0.694    187.256   
%     IRTNN   14.667    0.232    0.698    770.261   
%       TNN   12.989    0.285    0.681    21.971   
%     MF_TV   12.272    0.174    0.685    2210.062  
%    SNN_TV   14.742    0.559    0.671    122.428   
%    SPC_TV   15.765    0.421    0.715    441.207   
%    TNN_TV   18.087    0.648    0.681    87.273   
%      TCTV   26.584    0.767    0.869    121.428
%==========================================================================
%% Demo start
clc;
clear;close all;
rng('default');rng(1997);
addpath(genpath('lib'));
addpath(genpath('data'));
dataName = 'img_Lena';  % 'img_Einstein';
%  Please make sure the RGB image is a cubic of size [height, width, 3] and in range [0, 1].
%  You can use other tensor data such as Hyperspectral Image, Video, CT/MRI for test. 
%  Note some parameter might need reset for other methods. 
%  But the proposed TCTV is parameter free. 
dataRoad = ['data/', dataName];
resultsRoad = ['results/TC/results_for_', dataName];
if ~exist(resultsRoad); mkdir(resultsRoad); end
%% Set enable bits
Run_SNN       = 1;  % 1 or 0 <=> run or not
Run_BCPF      = 1;  % it occurs memory overflow for img "Einstein"
Run_KBR       = 1;
Run_IRTNN     = 1;
Run_TNN       = 1;
Run_MF_TV     = 1;
Run_SNN_TV    = 1;
Run_SPC_TV    = 1;
Run_TNN_TV    = 1;
Run_TCTV      = 1;  % our method

get_Recovered_Image = 1;  % if save the recovered image or not

%% Load Data 
methodName = {'Observed', 'SNN', 'BCPF', 'KBR', 'IRTNN', 'TNN', 'MF+TV', 'SNN+TV', 'SPC+TV', 'TNN+TV', 'TCTV'};
Mnum = length(methodName);
load(dataRoad);  % load data
[height, width, band] = size(data);
dim = [height, width, band];

%% Observation
i = 1;
missing_rate  = 0.9; % sampling rate, i.e, sampling_rate = 1 - missing rate

disp(['=== the missing rate is ', num2str(missing_rate), ' ===']);
saveRoad = ['results/TC/results_for_', dataName, '/', erase(num2str(missing_rate),'.')];
if ~exist(saveRoad); mkdir(saveRoad); end
if exist([saveRoad, '/QA_Results.mat']); load([saveRoad, '/QA_Results.mat']); end
if exist([saveRoad, '/Results.mat']); load([saveRoad, '/Results.mat']); end
if get_Recovered_Image
    if ~exist([saveRoad,'/Recovered_Image']); mkdir([saveRoad,'/Recovered_Image']); end
    imwrite(data, [saveRoad, '/Recovered_Image/Original.jpg']);
end

sampling_rate = 1-missing_rate;
m          = round(prod(dim)*sampling_rate);
sort_dim   = randperm(prod(dim));
Omega      = sort_dim(1:m); % sampling pixels' index
Obs        = zeros(dim);
Obs(Omega) = data(Omega); % observed Img

Results{i} = Obs;
[PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
enList = 1;

%% Run SNN
i = i+1;
if Run_SNN
    addpath(genpath(['competing methods/TC/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the HaLRTC.m or example.m in SNN file
    Omega_SNN = zeros(dim);
    Omega_SNN(Omega) = 1;
    Omega_SNN = logical(Omega_SNN);
    alpha = [1, 1, 1e-3];
    alpha = alpha / sum(alpha);
    maxIter = 500;
    epsilon = 1e-5;
    beta = 1e-2;

    tic;
    Results{i} = HaLRTC(data, Omega_SNN, alpha, beta, maxIter, epsilon);
    Time(i) = toc;
%     % or we can equally use the LRTC_SNN.m via ADMM sovler
%     addpath(genpath('competing methods/TC/other solvers'));
%     opts = [];
%     alpha = [15,15,1];
%     Results{i+1} = = LRTC_SNN(Obs,Omega,alpha,opts);
%     Time(i) = toc;
%     rmpath(genpath('competing methods/TC/other solvers'));
    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TC/', methodName{i}]));
    enList = [enList, i];
end

%% Run BCPF
i = i+1;
if Run_BCPF
    addpath(genpath(['competing methods/TC/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the BCPF_IC_denoise.m or DemoBayesCP.m in BCPF file
    Omega_BCPF = zeros(dim);
    Omega_BCPF(Omega) = 1;
    tic
    [model] = BCPF_IC(Obs, 'obs', Omega_BCPF, 'init', 'rand', 'maxRank', 100, 'maxiters', 20, ...
             'tol', 1e-4, 'dimRed', 1, 'verbose', 1);
    Results{i} = double(model.X);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TC/', methodName{i}]));
    enList = [enList, i];
end

%% Run KBR
i = i+1;
if Run_KBR
    addpath(genpath(['competing methods/TC/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the KBR_TC.m or Demo_TC_MSI.m in KBR file
    
    Omega_KBR     = zeros(dim);
    Omega_KBR(Omega) = 1;
    Omega_KBR     = (Omega_KBR > 0);
    
    Par_KBR.tol     = 1e-5;
    Par_KBR.maxIter = 1000;
    Par_KBR.maxSubiter = 1;%lambda
    Par_KBR.rho = 1.05;
    Par_KBR.mu = 1e-5;
    Par_KBR.lambda = 0.01; % 0.01 is even better
    tic
    Results{i} = KBR_TC(Obs, Omega_KBR, Par_KBR);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TC/', methodName{i}]));
    enList = [enList, i];
end

%% Run IRTNN
i = i+1;
if Run_IRTNN
    addpath(genpath(['competing methods/TC/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the IRTNN.m or run_IRTNN.m in KBR file
    IDX = Omega ;
    M = opRestriction(prod(dim), IDX);
    x = data(:);
    y = M(x,1);   
%         IRTNN's different penalties, here, we use the lp non-convex penalty
        fun = 'lp' ;           gamma = 0.5 ;
%         fun = 'scad' ;        gamma = 100 ;
%         fun = 'logarithm' ;   gamma = 10 ;
%         fun = 'mcp' ;         gamma = 10 ;
%         fun = 'cappedl1' ;    gamma = 1000 ;
%         fun = 'etp' ;         gamma = 0.1 ;
%         fun = 'geman' ;       gamma = 10 ;
%         fun = 'laplace' ;     gamma = 10 ;
    tic
    Results{i} = IRTNN(fun,y,M,height,width,band,gamma) ;
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TC/', methodName{i}]));
    enList = [enList, i];
end

%% Run TNN
i = i+1;
if Run_TNN
    addpath(genpath(['competing methods/TC/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the HTNN_FFT.m or Demo.m in TNN file
    % note that HTNN_FFT degrades into the TNN.m in Zhang Zemin's method in CVPR14/TSP16
    Omega_TNN     = zeros(dim);
    Omega_TNN(Omega) = 1;
    opts = [];
    opts.mu = 1e-4;
    opts.max_mu = 1e8;
    opts.max_iter =500;
    opts.DEBUG = 1;
    opts.rho = 1.2;
    opts.tol = 1e-6;
    tic
    Results{i} = HTNN_FFT(Obs, Omega_TNN, opts);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TC/', methodName{i}]));
    enList = [enList, i];
end

%% Run MF_TV
i = i+1;
if Run_MF_TV
    addpath(genpath(['competing methods/TC/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the LRTC_TV.m or Demo.m in MF_TV file
    known = Omega;
    Obs_data = data(known);
    [known, id]= sort(known);
    Obs_data= Obs_data(id);
    Y_tensor0= zeros(dim);
    Y_tensor0(known)= Obs_data;
    % Initialization of the factor matrices X and A
    ratio = 0.005;
    R=AdapN_Rank(data,ratio);
    for j = 1:3
        coNway(j) = prod(dim)/dim(j);
    end
    for j = 1:3
        Y0{j} = Unfold(Obs,dim,j);
        Y0{j} = Y0{j}';
        X0{j}= rand(coNway(j), R(j));
        A0{j}= rand(R(j),dim(j));
    end
    % run MF_TV
    opts=[];
    tic
    Results{i} = LRTC_TV(Y0, Obs_data, A0, X0, Obs, dim, known, opts, height, width);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TC/', methodName{i}]));
    enList = [enList, i];
end

%% Run SNN_TV
i = i+1;
if Run_SNN_TV
    addpath(genpath(['competing methods/TC/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the HTNN_FFT.m or Demo.m in TNN file
    index = zeros(m,3);
    for j=1:m
        [v1,v2,v3] = ind2sub(dim,Omega(j));
        index(j,:) = [v1,v2,v3];
    end
    value = Obs(Omega);
    lambda=0.02;
    N=3;
    alpha=[1/N, 1/N, 1/N];
    beta=[1,1,0];
    tsize=dim;
    tic
    Results{i} = LRTC_TV_I(index, value, lambda, alpha, beta, tsize, N );
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TC/', methodName{i}]));
    enList = [enList, i];
end

%% Run SPC_TV
i = i+1;
if Run_SPC_TV
    addpath(genpath(['competing methods/TC/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the SPC.m in SPC_TV file
    % note that HTNN_FFT degrades into the TNN.m in Zhang Zemin's method in CVPR14/TSP16
    T = Obs;
    Q = (T ~= 0);
    % parameters setting
    TVQV    = 'tv';        % 'tv' or 'qv' ;
    rho     = [0.01 0.01 0]; % smoothness (0.1 - 1.0) for 'qv' and (0.01 - 0.5) for 'tv' is recommended.
    % TVQV    = 'qv';        % 'tv' or 'qv' ;
    % rho     = [0.5 0.5 0]; % smoothness (0.1 - 1.0) for 'qv' and (0.01 - 0.5) for 'tv' is recommended.
    K       = 10;          % Number of components which are updated in one iteration. (typically 10)
    SNR     = 25;          % error bound
    nu      = 0.01;        % threshold for R <-- R + 1.
    maxiter = 1000;        % maximum number of iteration
    tol     = 1e-5;        % tolerance
    out_im  = 0;           % you can monitor the process of 'image' completion if out == 1. 'saved' directory is necessary to save the individual rank images.
    tic
    Results{i} = SPC(T,Q,TVQV,rho,K,SNR,nu,maxiter,tol,out_im);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TC/', methodName{i}]));
    enList = [enList, i];
end

%% Run TNN_TV
i = i+1;
if Run_TNN_TV
    addpath(genpath('competing methods/TC/other solvers'));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the LRTC_TNN_2DTV.m or Demo.m in 'other solvers' file
    opts = [];
    lambda = 10;
    tic
    Results{i} = LRTC_TNN_2DTV(Obs,Omega,lambda,opts);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    rmpath(genpath('competing methods/TC/other solvers'));
    enList = [enList, i];
end

%% Run TCTV
i = i+1;
if Run_TCTV
    addpath(genpath('TCTV'));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the TCTV_TC.m or Demo.m in TCTV file
    opts = [];
    opts.transform = 'DCT'; 
    % opts.transform = 'DFT'; % the DCT-based is usually better than the DFT-based
    opts.directions = [1 2]; % for images, [1 2]; for HSI, [1 2 3]
    tic
    Results{i} = TCTV_TC(Obs, Omega, opts);
    Time(i) = toc;
    [PSNR(i), SSIM(i), FSIM(i)] = Img_QA(data, Results{i});
    if get_Recovered_Image; imwrite(Results{i}, [saveRoad,'/Recovered_Image/', methodName{i}, '.jpg']); end
    addpath(genpath('TCTV'));
    enList = [enList, i];
end

%% Show result
fprintf('\n');    

fprintf('================== QA Results =====================\n');
fprintf(' %8.8s    %5.5s    %5.5s    %5.5s     %5.5s  \n',...
    'Method', 'MPSNR', 'MSSIM', 'MFSIM',  'Time');
% enList = 1:Mnum;
for i = 1:length(enList)
    fprintf(' %8.8s   %5.3f    %5.3f    %5.3f    %5.3f   \n',...
        methodName{enList(i)}, PSNR(enList(i)), SSIM(enList(i)), FSIM(enList(i)), Time(enList(i)));
end

fprintf('================== Show Results =====================\n');
figure;
plt_w = 5;
plt_h = ceil((length(enList)+1)/plt_w);
subplot(plt_h, plt_w,1); imshow(data);title('Original');
for ii = 1:length(enList)
    hold on;subplot(plt_h, plt_w,ii+1); imshow(Results{enList(ii)});title(methodName{enList(ii)});
end

%% save results
All = [PSNR; SSIM; FSIM; Time];
save([saveRoad,'/QA_Results'], 'All', 'PSNR', 'SSIM', 'FSIM', 'Time');
save([saveRoad,'/Results'], 'Results');



