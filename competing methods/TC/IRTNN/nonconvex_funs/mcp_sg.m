function y = mcp_sg(x,gamma,lambda)
%  supergradient of MCP
x = abs(x) ;
y = zeros(size(x)) ;
ind = find(x<gamma*lambda) ;
y(ind) = lambda-x(ind)/gamma ;