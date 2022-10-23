function y = etp_sg(x,gamma,lambda)
% supergradient of ETP
x = abs(x) ;
y = gamma*lambda*exp(-gamma*x)./(1-exp(-gamma)) ;
