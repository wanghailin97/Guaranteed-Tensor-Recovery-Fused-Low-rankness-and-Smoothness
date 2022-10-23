function y = laplace_sg(x,gamma,lambda)
% Laplace penalty

x = abs(x) ;
y = lambda*exp(-x/gamma)/gamma ;

