%==========================================================================
% This script compares multi-spectral imagery (MSI) completion methods
% listed as follows:
%   1. ADMM(ALM)-based matrix completion (MC-ALM)
%   2. HaLRTC
%   3. the factorization-based TC method (TMac)
%   4. joint trace/TV based TC method (Trace-TV)
%   5. tensor-SVD-based TC method (t-SVD)
%   6. minmax concave plus penalty-based TC method (McpTC)
%   7. smoothly clipped absolute deviation penalty-based TC method (ScadTC)
%   9. KBR-TC
%
% Four quality assessment (QA) indices -- PSNR, SSIM, FSIM, ERGAS
% -- are calculated for each methods after denoising.
%
% You can:
%       1. Type 'Demo_TC_MSI' to to run various methods and see the pre-computed results.
%       3. Use 'help Demo' for more information.
%       4. Change test MSI by simply modifying variable 'dataname' in Demo.m (NOTE: make sure your MSI
%          meets the format requirements).
%       5. Change sampling rate by modifying variables  'sample_ratio ' in Demo_TC_MSI.m
%       6. Select competing methods by turn on/off the enable-bits in Demo_TC_MSI.m
%
% See also KBR_TC, Tc_liu, Mc_adm, TMac, LRTV_TC, McpLRTC, SCADLRTC, tensor_cpl_admm and MSIQA
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
dataname = 'testMSI_2';%  Please make sure this MSI is of size height x width x nbands and in range [0, 1].
                       %  can also use 'testMSI_1' and 'testMSI_1', as other examples.
dataRoad = ['data/' dataname];
saveroad = ['result/result_for_' dataname];

%% Set enable bits
sample_ratio = 0.1; % sampling ratio;  % higher sigma_ratio <--> heavier Gaussian noise
disp( '=== The variance of noise is 0.1 ===');
EN_MC_ALM   = 1;  % set to 0 for turning off;
EN_HaLRTC   = 1;
EN_Tmac     = 1;
EN_Trace_TV = 1;
EN_t_SVD    = 1;
EN_McpTC    = 1;
EN_ScadTC   = 1;
EN_KBR_TC   = 1;
getGif      = 1; % whether save gif result or not
getImage    = 1; % whether save the appointed band of the reconstructed MSI as image result or not
mkdir(saveroad);
rng(1);
%% initial Data
methodname  = {'MC-ALM', 'HaLRTC',  'Tmac', 'Trace-TV','t-SVD', 'McpTC', 'ScadTC','KBR-TC'};
Mnum   = length(methodname);
load(dataRoad); % load data
T            = normalized(Omsi);
band2show    = 15; %the band to show and save
temp         = T(:,:,band2show);
maxI         = max(temp(:));
minI         = min(temp(:));
sizeT        = size(T);
ndim         = length(sizeT);

%% initialization of the parameters
alpha   = ones(1, 3);
alpha   = alpha / sum(alpha);
maxIter = 1000;
epsilon = 1e-5;
beta    = 1e-2; %for the ||_{*}||_{*}||_{*}

%%  sampling with random position
Omega     = zeros(sizeT);
ind       = randperm(prod(sizeT));
known     = ind(1:floor(prod(sizeT)*sample_ratio));
Omega(known) = 1;
Omega     = (Omega > 0);
T_miss = T;
T_miss(~Omega) = mean(T(Omega));
enList = [];

if getGif; mkdir([saveroad, '/GIF']); end;
if getImage; mkdir([saveroad, '/Image']);end


%% Use MC_ALM
i  = 1;
if EN_MC_ALM
    disp(['performing ',methodname{i}, ' ... ']);
    tic;
    Re_msi{i} = Mc_adm(T, Omega, beta, maxIter, epsilon);
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(T * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band2show)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end



%% Use HaLRTC
i = i+1;
if EN_HaLRTC
    disp(['performing ',methodname{i}, ' ... ']);
    tic;
    [Re_msi{i}] = Tc_liu(...
        T, ...                       % a tensor whose elements in Omega are used for estimating missing value
        Omega,...                    % the index set indicating the obeserved elements
        alpha,...                    % the coefficient of the objective function, i.e., \|X\|_* := \alpha_i \|X_{i(i)}\|_*
        beta,...                     % the initial value of the parameter; it should be small enough
        maxIter,...                  % the maximum iterations
        epsilon...                   % the tolerance of the relative difference of outputs of two neighbor iterations
        );
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(T * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band2show)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end

%% Use Tmac
i = i+1;
if EN_Tmac
    disp(['performing ',methodname{i}, ' ... ']);
    opts = [];
    opts.maxit = maxIter;
    opts.tol = epsilon; % run to maxit by using negative tolerance
    opts.Mtr = T; % pass the true tensor to calculate the fitting
    opts.alpha_adj = 1;
    opts.rank_adj = -ones(3,1);
    opts.rank_min = [160,160,1];
    opts.rank_max = [160,160,6];
    EstCoreNway   = round(([0.7*sizeT(1), 0.7*sizeT(2), 1]+[sizeT(1), sizeT(2), 10])/2);
    coNway = zeros(1,ndim);
    for n = 1:ndim
        coNway(n) = prod(sizeT)/sizeT(n);
    end
    % use random generated starting point
    for k = 1:ndim
        U0{k} = randn(sizeT(k),EstCoreNway(k));
        V0{k} = randn(EstCoreNway(k),coNway(k));
    end
    %         known   = find(Omega>eps);
    opts.X0 = U0; opts.Y0 = V0;
    tic;
    [U_dec,V_dec,Out_dec] = TMac(T_miss(known),known,sizeT,EstCoreNway,opts);
    Re_msi{i} = zeros(sizeT);
    for k = 1:ndim
        Re_msi{i} = Re_msi{i}+Out_dec.alpha(k)*Fold(U_dec{k}*V_dec{k},sizeT,k);
    end
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(T * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band2show)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end

%% Use Trace/TV
i = i+1;
if EN_Trace_TV
    disp(['performing ',methodname{i}, ' (will take about 100s)... ']);
    tic;
    Re_msi{i} = LRTV_TC(T,Omega);
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(T * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band2show)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end


%% Use t_SVD
i = i+1;
if EN_t_SVD
    disp(['performing ',methodname{i}, ' (will take about 400s)... ']);
    tic;
     normalize              =        max(T_miss(:))                ;
        inX                    =        T_miss/normalize              ;
        gamma                  =        1                             ;
        maxItr                 =        maxIter                       ; % maximum iteration
        rho                    =        0.01                          ;
        myNorm                 =        'tSVD_1'                      ; % dont change for now
        A                      =        diag(sparse(double(Omega(:)))); % sampling operator
        b                      =        A * inX(:)                    ; % available data
        Re_msi{i}    =    tensor_cpl_admm( A , b , rho , gamma , ...
            sizeT , maxItr , myNorm , 1 );
        Re_msi{i}                   =        Re_msi{i} * normalize    ;
        Re_msi{i}                   =        reshape(Re_msi{i},sizeT) ;     
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(T * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band2show)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end

%% Use McpTC
i = i+1;
if EN_McpTC
    disp(['performing ',methodname{i}, ' (will take about 600s)... ']);
    tic;
    Re_msi{i}  = McpLRTC(T,Omega,alpha,beta,maxIter,epsilon);
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(T * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band2show)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end

%% Use ScadTC
i = i+1;
if EN_ScadTC
    disp(['performing ',methodname{i}, ' (will take about 600s)... ']);
    tic;
    Re_msi{i} = SCADLRTC(T,Omega,alpha,beta,maxIter,epsilon);
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(T * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band2show)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end

%% Use KBR_TC
i = i+1;
if EN_KBR_TC
    disp(['performing ',methodname{i}, ' (will take about 400s)... ']);
    Par2.tol     = epsilon;
    Par2.maxIter = maxIter;
    Par2.maxSubiter = 1;%lambda
    Par2.rho = 1.05;
    Par2.mu = beta*1e-3;
    Par2.lambda = 0.1; % 0.01 is even better
    tic;
    Re_msi{i} = KBR_TC(T.*Omega, Omega, Par2);
    Time(i) = toc;
    [psnr(i), ssim(i), fsim(i), ergas(i)] = MSIQA(T * 255, Re_msi{i}  * 255);
    disp([methodname{i}, ' done in ' num2str(Time(i)), ' s.'])
    disp([methodname{i}, ' PSNR ' num2str(psnr(i)), ' .'])
    disp('...')
    if getGif; togetGif(Re_msi{i},[saveroad, '/GIF/', methodname{i}]); end;
    if getImage; imwrite((Re_msi{i}(:,:,band2show)-minI)/(maxI-minI),[saveroad,'/Image/', methodname{i}, '.png']);end
    enList = [enList,i];
end

%% Show result
fprintf('\n');
fprintf('================== Result =====================\n');
fprintf(' %8.8s    %5.4s    %5.4s    %5.4s    %5.5s    \n','method','PSNR', 'SSIM', 'FSIM', ' ERGAS');
for i = 1:length(enList)
    fprintf(' %8.8s    %5.3f    %5.3f    %5.3f    %5.3f    \n',...
        methodname{enList(i)},psnr(enList(i)), ssim(enList(i)), fsim(enList(i)), ergas(enList(i)));
end
fprintf('================== Result =====================\n');
showMSIResult(Re_msi,T,T_miss, methodname,enList,band2show);
save([saveroad,'\Result'], 'psnr','ssim','fsim','ergas','methodname','Re_msi');