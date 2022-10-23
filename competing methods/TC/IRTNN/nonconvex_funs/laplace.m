function y = laplace(x,gamma,lambda)
% Laplace penalty

x = abs(x) ;
y = lambda*(1-exp(-x/gamma)) ;

