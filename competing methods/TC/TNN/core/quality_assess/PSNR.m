function p = PSNR(x,y)

% psnr - compute the Peack Signal to Noise Ratio, defined by :
%       PSNR(x,y) = 10*log10( max(max(x),max(y))^2 / |x-y|^2 ).
%
%   p = psnr(x,y);
%
%   Copyright (c) 2004 Gabriel Peyr

d = mean( mean( (x(:)-y(:)).^2 ) );
m1 = max( abs(x(:)) );
m2 = max( abs(y(:)) );
% m = max(m1,m2);
m = 1;

p = 10*log10( m^2/d );

%  p = 10*log10(1/mean((x(:)-y(:)).^2));
% p=10*( log10( 255^2 ) + (log10(size(x,1)*size(x,2))-log10(norm(x(:)-y(:))^2)));
