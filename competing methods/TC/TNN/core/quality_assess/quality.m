function [psnr, ssim, fsim] = quality(imagery1, imagery2)
%==========================================================================
%
% Input:
%   imagery1 - the reference tensor
%   imagery2 - the target tensor

%
% Output:
%   psnr - Peak Signal-to-Noise Ratio
%   ssim - Structure SIMilarity
%   fsim - Feature SIMilarity
%==========================================================================


Nway = size(imagery1);
psnr = zeros(prod(Nway(3:end)),1);

ssim = psnr;
fsim = psnr;
for i = 1:prod(Nway(3:end))   
    psnr(i) = psnr_index(imagery1(:, :, i), imagery2(:, :, i));
    ssim(i) = ssim_index(imagery1(:, :, i), imagery2(:, :, i));
    fsim(i) = FeatureSIM(imagery1(:, :, i), imagery2(:, :, i));
end
psnr = nanmean(psnr);
ssim = nanmean(ssim);
fsim = nanmean(fsim);

