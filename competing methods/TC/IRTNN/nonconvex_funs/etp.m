function y = etp(x,gamma,lambda)
% ETP
x = abs(x) ;
y = lambda*(1-exp(-gamma*x))./(1-exp(-gamma)) ;
