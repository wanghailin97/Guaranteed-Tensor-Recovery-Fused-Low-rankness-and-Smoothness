function y = cappedl1_sg(x,gamma,lambda)
% supergradient of Capped L1 penalty
x = abs(x) ;
y = zeros(size(x)) ;
ind = x<gamma ;
y(ind) = lambda ;
