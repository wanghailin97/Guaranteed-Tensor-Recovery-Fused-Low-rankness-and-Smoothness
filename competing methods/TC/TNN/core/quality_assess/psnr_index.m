function p = psnr_index(x,y)

% psnr - compute the Peack Signal to Noise Ratio, defined by :
p=10*log10(255^2/mse(x-y));



