function y = lp(x,p,lambda)
% lp penalty

x = abs(x) ;
epsilon = 0 ;
y = lambda*(x+epsilon).^p ;

