function y = geman_sg(x,gamma,lambda)
% Geman penalty

x = abs(x) ;
y = lambda*gamma./(x+gamma).^2 ;

