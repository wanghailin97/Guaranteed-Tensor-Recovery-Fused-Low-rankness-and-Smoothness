%==========================================================================
% This script compares background subtraction  methods
% listed as follows:
%   1. ADMM(ALM)-based RPCA (IALM-RPCA)
%   2. HoRPCA
%   3. tensor-SVD-based RPCA method (t-SVD)
%   4. ITS-RPCA
%
% Four quality assessment (QA) indices -- Fmeasure
% -- are calculated for each methods after denoising.
%
% You can:
%       1. Type 'Demo_RPCA_video' to to run various methods and see the pre-computed results.
%       2. Change test MSI by simply modifying variable 'dataname' in Demo_RPCA_video.m (NOTE: make sure your MSI
%          meets the format requirements).
%       3. Select competing methods by turn on/off the enable-bits in Demo_RPCA_video.m
%
% See also KBR_RPCA, inexact_alm_rpca, tensor_rpca_adal, tensor_rpca  and findFMeasure
%
% more detail can be found in [1]
% [1] Qi Xie, Qian Zhao, Deyu Meng*, Zongben Xu. 
%     Kronecker-Basis-Representation Based Tensor Sparsity and Its Applications to Tensor Recovery,  
%     IEEE Transactions on Pattern Analysis and Machine Intelligence, 2017.
%
% by Qi Xie 
% 10/24/2017
%==========================================================================
clc;
clear;close all;
addpath(genpath('lib'));
dataname = 'testVideo';%  Please make sure this MSI is of size height x width x nbands and in range [0, 1].
truthname = 'groundtruth_at_16th_frame.bmp';
dataRoad = ['data/' dataname];
saveroad = ['result/result_for_' dataname];

%% Set enable bits
disp( '=== The variance of noise is 0.1 ===');
EN_IALM_RPCA   = 1;  % set to 0 for turning off;
EN_HoRPCA      = 1;
EN_t_SVD       = 1;
EN_ITS_RPCA    = 1;
getGif         = 0; % whether save gif result or not
mkdir(saveroad);
rng(1);
%% initial Data
methodname  = {'IALM-RPCA','HoRPCA','t-SVD','KBR-RPCA'};
Mnum   = length(methodname);
load(dataRoad); % load data
S0 = imread(['data/',truthname]);
S0 = (double(mean(S0,3))>100);
testNum = 16;  % the frame number of the groundtruth frame

%% initialization of the parameters
alpha   = ones(1, 3);
alpha   = alpha / sum(alpha);
maxIter = 1000;
mu      = 1e-2; 
epsilon = 1e-5;

sizeD     = size(D);
ndim      = length(sizeD);

%% normalization
normalization = 255;
D             = D/normalization;

enList = [];

%% Use IALM-RPCA
i  = 1;
if EN_IALM_RPCA
    disp(['performing ',methodname{i}, ' ... ']);
    tic;
    [B_ALm, T_ALm] = inexact_alm_rpca(Unfold(D,sizeD,3)',0.4 / sqrt(sizeD(2)*sizeD(1)));
    B{i}       = Fold(B_ALm',sizeD,3);
    F{i}       = Fold(T_ALm',sizeD,3);
    Time(i) = toc;
    [fmeasure(i), S{i}, bestTH(i)] = findFMeasure(F{i}(:,:,testNum), S0);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' Fmeasure ' num2str(fmeasure(i)), ' .'])
    disp('...')
    if getGif; togetGif(B{i},[saveroad, '/GIF_Background_of_/', methodname{i}]); end;
    if getGif; togetGif(F{i},[saveroad, '/GIF_Foreground_of_/', methodname{i}]); end;
    enList = [enList,i];
end

%% Use HoRPCA
i  = i + 1;
if EN_HoRPCA
    disp(['performing ',methodname{i}, ' ... ']);
    
    ParH.lambda      = .3/sqrt(max(sizeD));
    ParH.mu1         = 5*std(D(:));
    ParH.mu2         = 5*std(D(:));
    ParH.max_iter    = 300;
    ParH.verbose     = false;
    data.T           = D;
    data.X           = D;
    ParH.X0          = B{1};
    ParH.E0          = F{1};
    ParH.V0          = cell(3,1);
    ParH.opt_tol     = 1e-5;
    for vi = 1:3
        ParH.V0{vi}  = zeros(size(F{1}));
    end
    tic;
    resultsHOrpca    = tensor_rpca_adal(data, ParH);
    B{i}       = double(resultsHOrpca.X);
    F{i}       = double(resultsHOrpca.E);
    Time(i) = toc;
    [fmeasure(i), S{i}, bestTH(i)] = findFMeasure(F{i}(:,:,testNum), S0);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' Fmeasure ' num2str(fmeasure(i)), ' .'])
    disp('...')
    if getGif; togetGif(B{i},[saveroad, '/GIF_Background_of_/', methodname{i}]); end;
    if getGif; togetGif(F{i},[saveroad, '/GIF_Foreground_of_/', methodname{i}]); end;
    enList = [enList,i];
end

%% Use t-SVD
i  = i + 1;
if EN_t_SVD
    disp(['performing ',methodname{i}, ' ... ']);
    tic;
    [B{i}, F{i}] =   tensor_rpca( D , .03/sqrt(size(D,1)));
    Time(i) = toc;
    [fmeasure(i), S{i}, bestTH(i)] = findFMeasure(F{i}(:,:,testNum), S0);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' Fmeasure ' num2str(fmeasure(i)), ' .'])
    disp('...')
    if getGif; togetGif(B{i},[saveroad, '/GIF_Background_of_/', methodname{i}]); end;
    if getGif; togetGif(F{i},[saveroad, '/GIF_Foreground_of_/', methodname{i}]); end;
    enList = [enList,i];
end

%% Use KBR-RPCA
i  = i + 1;
if EN_ITS_RPCA
    beta           = 2.5*sqrt(max(sizeD));
    gamma          = beta*100;
    Par.maxIter    = 1000;
    Par.lambda     = 0.1;
    Par.mu         = mu*1000;
    Par.tol        = 1e-5;
    Par.rho        = 1.05;
    disp(['performing ',methodname{i}, ' ... ']);
    tic;
    [B{i}, F{i}] =   KBR_RPCA(D,beta,gamma,Par);
    Time(i) = toc;
    [fmeasure(i), S{i}, bestTH(i)] = findFMeasure(F{i}(:,:,testNum), S0);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' Fmeasure ' num2str(fmeasure(i)), ' .'])
    disp('...')
    if getGif; togetGif(B{i},[saveroad, '/GIF_Background_of_/', methodname{i}]); end;
    if getGif; togetGif(F{i},[saveroad, '/GIF_Foreground_of_/', methodname{i}]); end;
    enList = [enList,i];
end
save([saveroad,'\Result'], 'fmeasure','B','S','F');
%% Show result
fprintf('\n');
fprintf('================== Result =====================\n');
fprintf(' %9.9s    %8.8s    \n','method','Fmeasure');
for i = 1:length(enList)
    fprintf(' %9.9s    %4.2f     \n',...
        methodname{enList(i)},fmeasure(enList(i)));
end
fprintf('================== Result =====================\n');
showVideoResult(D,B, methodname,enList,35)