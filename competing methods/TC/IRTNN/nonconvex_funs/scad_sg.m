function y = scad_sg(x,gamma,lambda)
%  supergradient of SCAD

x = abs(x) ;

y = zeros(size(x)) ;
ind = x<=lambda ;
y(ind) = lambda ;
ind = find(x>lambda & x<=gamma*lambda) ;
y(ind) = (gamma*lambda-x(ind))/(gamma-1) ;

