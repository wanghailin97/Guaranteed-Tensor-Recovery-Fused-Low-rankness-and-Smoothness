function y = cappedl1(x,gamma,lambda)
% Capped L1 penalty
x = abs(x) ;
y = ones(size(x))*gamma*lambda ;
ind = find(x<gamma) ;
y(ind) = lambda*x(ind) ;
