function y = logarithm_sg(x,gamma,lambda)
% logarithm penalty

x = abs(x) ;
y = lambda*gamma./((gamma*x+1).*log(gamma+1)) ;