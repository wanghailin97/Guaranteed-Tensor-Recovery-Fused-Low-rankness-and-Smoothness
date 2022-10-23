function y = geman(x,gamma,lambda)
% Geman penalty

x = abs(x) ;
y = lambda*x./(x+gamma) ;

