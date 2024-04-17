%% Demo of Tensor Robust Principle Component Analysis
%  test for TCTV based Hyperspectral image denoising and comparison with other SOTA methods
%==========================================================================
% This script compares many tensor completion methods
% listed as follows:
%   1. SNN
%   2. KBR
%   3. TNN
%   4. LRTV
%   5. LRTDTV
%   6. TLR_HTV
%   7. TCTV (proposed method)
%
% -- five quality assessment (QA) indices -- MPSNR, MSSIM, MFSIM, RRGAS, MSAM
% -- are calculated for each method's denoising result, see HSI_QA
%
% You can:
%       1. Type 'Demo' to to run various methods and see the pre-computed results.
%       2. Change test HSI by simply modifying variable 'dataName' in Demo.m 
%          (NOTE: make sure your HSI meets the format requirements).
%       3. Change sparse corruption level by modifying variables  'sparse_rate' in Demo.m
%       4. Select competing methods by turn on/off the enable-bits in Demo.m
%
% More detail can be found in 
% Hailin Wang et.al, Guaranteed Tensor Recovery Fused Low-rankness and Smoothness, 2022
%
% by Hailin Wang, 2022
%==========================================================================
% The denoising results of the HSI "PaC" when sparse_rate = 20%
% ================== QA Results =====================
%    Method    MPSNR    MSSIM    MFSIM    ERGAS     MSAM     Time  
%     Noisy   11.512    0.128    0.578    989.25    45.48    0.000   
%       SNN   41.136    0.995    0.995    45.726    1.377    183.52   
%       KBR   32.774    0.961    0.976    90.501    4.973    71.257   
%       TNN   44.477    0.987    0.991    49.203    4.175    80.186   
%      LRTV   36.695    0.977    0.987    55.226    3.038    21.910   
%    LRTDTV   37.702    0.978    0.990    48.689    3.270    69.091   
%   TLR_HTV   34.816    0.963    0.977    71.906    4.075    32.791   
%      TCTV   46.461    0.991    0.994    36.614    3.402    352.36 
%==========================================================================
% The denoising results of the HSI "PaC" when sparse_rate = 40%
% ================== QA Results =====================
%    Method    MPSNR    MSSIM    MFSIM    ERGAS     MSAM     Time  
%     Noisy   8.501     0.048    0.434    1398.3    48.76    0.000   
%       SNN   36.221    0.983    0.988    75.095    2.140    181.80   
%       KBR   30.497    0.939    0.958    112.64    5.942    69.683   
%       TNN   36.835    0.966    0.979    69.883    6.310    81.770   
%      LRTV   33.156    0.949    0.966    81.819    4.081    20.821   
%    LRTDTV   35.671    0.968    0.983    61.667    3.909    66.364   
%   TLR_HTV   31.447    0.923    0.950    101.58    5.397    31.596   
%      TCTV   41.890    0.986    0.991    45.287    4.008    155.01  
%==========================================================================
%% 
clc;
clear;close all;
rng('default');rng(1997);
addpath(genpath('lib'));
addpath(genpath('data'));
dataName = 'hsi_PaC'; 
%  Please make sure the HSI is a cubic of size [height, width, band] and in range [0, 1].
%  You can use other tensor data such as RGB Image, Video, CT/MRI for test. 
%  Note some parameter might need reset for other methods. 
%  But the proposed TCTV is parameter free, where the trade-off paramter 'lambda' is determined theoretically 
dataRoad = ['data/', dataName];
resultsRoad = ['results/TRPCA/results_for_', dataName];
if ~exist(resultsRoad); mkdir(resultsRoad); end
%% Set enable bits
Run_SNN       = 1;  % 1 or 0 <=> run or not
Run_KBR       = 1;
Run_TNN       = 1;
Run_LRTV      = 1;
Run_LRTDTV    = 1;
Run_TLR_HTV   = 1;
Run_TCTV      = 1;  % our method

getGIF         = 1;  % if save gif result or not
getBand        = 1;  % if save one slected band of the denoised HSI as a gray image result or not
selected_band  = 25;
getPseudoImage = 1;  % if save three slected bands of the denoised HSI as a pseudo color image result or not
selected_bands = [49, 27, 7];

%% Load Data 
methodName = {'Noisy', 'SNN',  'KBR', 'TNN', 'LRTV', 'LRTDTV', 'TLR-HTV','TCTV'};
Mnum = length(methodName);
load(dataRoad);  % load data
Ohsi = data;
[height, width, band] = size(Ohsi);
dim = [height, width, band];

%% Add impluse noise
i = 1;
sparse_rate  = 0.4; % the ratio of salt&pepper noise 

disp(['=== the noise level is ', num2str(sparse_rate), ' ===']);
saveRoad = ['results/TRPCA/results_for_', dataName, '/', erase(num2str(sparse_rate),'.')];
if ~exist(saveRoad); mkdir(saveRoad); end
if exist([saveRoad, '/QA_Results.mat']); load([saveRoad, '/QA_Results.mat']); end
if exist([saveRoad, '/Results.mat']); load([saveRoad, '/Results.mat']); end
if getGIF
    if ~exist([saveRoad, '/GIF']); mkdir([saveRoad, '/GIF']); end
    togetGif(Ohsi, [saveRoad, '/GIF/Ohsi']); 
end
if getBand
    if ~exist([saveRoad,'/Band']); mkdir([saveRoad,'/Band']); end
    imwrite(Ohsi(:,:,selected_band), [saveRoad, '/Band/Ohsi.jpg']); 
end
if getPseudoImage
    if ~exist([saveRoad,'/PseudoImage']); mkdir([saveRoad,'/PseudoImage']); end
    imwrite(PseudoImage(Ohsi, selected_bands), [saveRoad, '/PseudoImage/Ohsi.jpg']);
end

for ii = 1:band
    Nhsi(:,:,ii) = imnoise(Ohsi(:,:,ii), 'salt & pepper', sparse_rate);
end

Results{i} = Nhsi;
[MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
enList = 1;

%% Run SNN
i = i+1;
if Run_SNN
    addpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']); 
    % see the trpca_snn.m in SNN file
    % the used code is from Lu Canyi, A unified alternating direction method of multipliers by majorization minimization, TPAMI, 2018 
    opts = [];
    opts.DEBUG = 1;
    tic
    alpha=[1, 1, 200]; 
    Results{i} = trpca_snn(Nhsi, alpha, opts);
    Time(i) = toc;
    [MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
    if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
    if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
    if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    enList = [enList, i];
end

%% Run KBR
i = i+1;
if Run_KBR
    addpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the KBR_RPCA.m or Demo_TC_MSI.m in KBR file
    
    beta           = 2.5*sqrt(max(dim));
    gamma          = beta*100;
    Par.maxIter    = 1000;
    Par.lambda     = 0.1;
    Par.mu         = 10;
    Par.tol        = 1e-5;
    Par.rho        = 1.1;
    
    tic
    Results{i} =   KBR_RPCA(Nhsi,beta,gamma,Par);
    Time(i) = toc;
    [MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
    if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
    if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
    if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    enList = [enList, i];
end

%% Run TNN
i = i+1;
if Run_TNN
    addpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the trpca_tnn.m  in TNN file
    lambda = 1/sqrt(min(height,width)*band);
    %lambda = 0.05;
    opts = [];
    opts.DEBUG = 1;
    
    tic
    Results{i} = trpca_tnn(Nhsi, lambda, opts);
    Time(i) = toc;
    [MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
    if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
    if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
    if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    enList = [enList, i];
end

%% Run LRTV
i = i+1;
if Run_LRTV
    addpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the LRTV_accelerate.m in LRTV file
    tau = 0.005;
    lambda = 5/sqrt(height*width);
    rank = 8;
    
    tic
    Results{i} = LRTV_accelerate(Nhsi, tau, lambda, rank);
    Time(i) = toc;
    [MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
    if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
    if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
    if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    enList = [enList, i];
end

%% Run LRTDTV
i = i+1;
if Run_LRTDTV
    addpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the LRTDTV.m or in TRTDTV file
    tau    = 1;
    lambda = 10;
    Rank   = [120, 120, 8];
    
    tic
    Results{i} = LRTDTV(Nhsi, tau, lambda, Rank);
    Time(i) = toc;
    [MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
    if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
    if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
    if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    enList = [enList, i];
end

%% Run TLR_HTV
i = i+1;
if Run_TLR_HTV
    addpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the SPC.m in SPC_TV file
    % note that HTNN_FFT degrades into the TNN.m in Zhang Zemin's method in CVPR14/TSP16
    Omega = ones(dim);
    
    opts = [];
    opts.theta  = 100;
    opts.xi     = 1.2;
    opts.tau    = 1.618;
    opts.rho    = 1000;
    opts.lambda = 2/sqrt(max(height,width)*band);
    %opts.lambda = 0.05;
    opts.beta   = 0.015; %rho
    opts.gamma  = 0.005;
    opts.tol    = 1e-3;
    opts.itemax = 30;
    opts.DEBUG = 1;
    tic
    Results{i} = LowSparN(Nhsi, Omega, opts);
    Time(i) = toc;
    [MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
    if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
    if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
    if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
    rmpath(genpath(['competing methods/TRPCA/', methodName{i}]));
    enList = [enList, i];
end

%% Run TCTV
i = i+1;
if Run_TCTV
    addpath(genpath('TCTV'));
    disp(['Running ',methodName{i}, ' ... ']);
    % use the recommended parameters from its released code 
    % see the TCTV_TRPCA.m in TCTV file
    opts = [];
    opts.rho = 1.25; % larger rho makes the algorithm faster while lose the accuracy
    opts.directions = [1,2,3]; % consider the lowrankness and smoothness both along the spatial and spectral directions
    tic
    Results{i} = TCTV_TRPCA(Nhsi, opts);
    Time(i) = toc;
    [MPSNR(i), MSSIM(i), MFSIM(i), ERGAS(i), MSAM(i)] = HSI_QA(Ohsi * 255, Results{i} * 255);
    if getGIF; togetGif(Results{i}, [saveRoad, '/GIF/', methodName{i}]); end
    if getBand; imwrite(Results{i}(:,:,selected_band), [saveRoad,'/Band/', methodName{i}, '.jpg']); end
    if getPseudoImage; imwrite(PseudoImage(Results{i}, selected_bands), [saveRoad,'/PseudoImage/', methodName{i}, '.jpg']); end
    rmpath(genpath('TCTV'));
    enList = [enList, i];
end

%% Show result
fprintf('\n');

fprintf('================== QA Results =====================\n');
fprintf(' %8.8s    %5.5s    %5.5s    %5.5s    %5.5s    %5.5s    %5.5s  \n',...
    'Method', 'MPSNR', 'MSSIM', 'MFSIM', 'ERGAS', 'MSAM', 'Time');
for i = 1:length(enList)
    fprintf(' %8.8s   %5.3f    %5.3f    %5.3f    %5.3f    %5.3f    %5.3f   \n',...
        methodName{enList(i)}, MPSNR(enList(i)), MSSIM(enList(i)), MFSIM(enList(i)), ERGAS(enList(i)), MSAM(enList(i)), Time(enList(i)));
end

fprintf('================== Show Results =====================\n');
close all;
showHSIResult(Results,Ohsi,methodName,enList,selected_band,band);


%% save results
All = [MPSNR; MSSIM; MFSIM; ERGAS; MSAM; Time];
save([saveRoad,'/QAResults'], 'All', 'MPSNR', 'MSSIM', 'MFSIM', 'ERGAS', 'MSAM', 'Time');
save([saveRoad,'/Results'], 'Results');



