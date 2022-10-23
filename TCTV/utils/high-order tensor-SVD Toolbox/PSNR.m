function psnr = PSNR(Xfull,Xrecover,maxP)
% 
% Written by  Wenjin Qin  (qinwenjin2021@163.com)
%

Xrecover = max(0,Xrecover);
Xrecover = min(maxP,Xrecover);

Xfull = max(0,Xfull);
Xfull = min(maxP,Xfull);

MSE = (norm(Xfull(:)-Xrecover(:))^2)/ numel(Xrecover);
psnr = 10*log10(maxP^2/MSE);