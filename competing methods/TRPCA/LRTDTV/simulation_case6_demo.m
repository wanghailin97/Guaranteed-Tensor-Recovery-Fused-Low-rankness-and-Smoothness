clear all;
clc;	

addpath ./prox_operators
addpath ./mylib
addpath ./quality_assess

% simulated experiment 6
%% load image
load simu_indian
Omsi = simu_indian;
load Simu_ratio                       
load Simu_noiselevel                 

Nmsi      = Omsi;
[M,N,p]   = size(Omsi);
%% Gaussian noise
for i = 1:p
     Nmsi(:,:,i)=Omsi(:,:,i)  + noiselevel(i)*randn(M,N);
end
%% S&P noise
for i = 1:p
     Nmsi(:,:,i)=imnoise(Nmsi(:,:,i),'salt & pepper',ratio(i));
end
%% dead lines 
for i=91:130
    indp=randperm(10,1)+2;
    ind=randperm(N-1,indp);
    an=funrand(2,length(ind));
    % searching the location of an which value is 1,2,3
    loc1=find(an==1);loc2=find(an==2);loc3=find(an==3);
    Nmsi(:,ind(loc1),i)=0; 
    Nmsi(:,ind(loc2):ind(loc2)+1,i)=0;
    Nmsi(:,ind(loc3)-1:ind(loc3)+1,i)=0;
end  
%% stripe line
for band=161:190
    num = 19+randperm(21,1);
    loc = ceil(N*rand(1,num));
    t = rand(1,length(loc))*0.5-0.25;
    Nmsi(:,loc,band) = bsxfun(@minus,Nmsi(:,loc,band),t);
end
%% LRTDTV denoising
tau    = 1;
lambda = 10;
Rank   = [120,120,10];
it     = 1;
clean_image                     = LRTDTV(Nmsi, tau,lambda,Rank);
[mpsnr(it),mssim(it),ergas(it)] = msqia(Omsi, clean_image);
it=2;
weight = [1,1,0.6];
clean_image                    = LRTDTV_w(Nmsi, tau,lambda,Rank,weight);
[mpsnr(it),mssim(it),ergas(it)]= msqia(Omsi, clean_image);