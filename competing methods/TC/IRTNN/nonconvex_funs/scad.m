function y = scad(x,gamma,lambda)
% SCAD

x = abs(x) ;
y = ones(size(x))*(gamma+1)*lambda^2/2 ;
ind = find(x<=lambda) ;
y(ind) = lambda*x(ind) ;
ind = find(x>lambda & x<=gamma*lambda) ;
y(ind) = (-x(ind).^2+2*gamma*lambda*x(ind)-lambda^2)./(2*(gamma-1)) ;

