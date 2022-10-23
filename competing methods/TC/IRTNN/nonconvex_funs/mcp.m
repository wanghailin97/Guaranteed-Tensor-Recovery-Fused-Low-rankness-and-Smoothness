function y = mcp(x,gamma,lambda)
% MCP
x = abs(x) ;
y = ones(size(x))*gamma*lambda^2*0.5 ;
ind = find(x<gamma*lambda) ;
y(ind) = lambda*x(ind)-x(ind).^2/(2*gamma) ;