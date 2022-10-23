function y = logarithm(x,gamma,lambda)
% logarithm penalty

x = abs(x) ;
y = lambda*log(gamma*x+1)./log(gamma+1) ;
