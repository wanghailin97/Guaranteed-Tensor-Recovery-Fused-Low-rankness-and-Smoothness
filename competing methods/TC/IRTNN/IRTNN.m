function X = IRTNN(fun,y,M,n1,n2,n3,gamma,err,x_initial,normfac,insweep,tol,decfac)

hfun_sg = str2func( [ fun '_sg'] ) ;

if nargin < 7
    gamma = 1;
end
if nargin < 8
    err = 1e-6;
end
if nargin < 9
    x_initial = zeros(n1*n2*n3,1);
end
if nargin < 10
    normfac = 1;
end
if nargin < 11
    insweep = 200;
end
if nargin < 12
    tol = 1e-5;    
end
if nargin < 13
    decfac = 0.9;
end

mu = 1.1*normfac;
x = x_initial;
lambdaInit = decfac*max(abs(M(y,2))); lambda = lambdaInit;
f_current = norm(y-M(x,1)) + lambda*norm(x,1);
while lambda > lambdaInit*tol
%     obj = [] ;
    for ins = 1:insweep    
%         obj(ins) = sum(hfun(svd(reshape(x,sizeX),'econ'),gamma,lambda)) + 0.5*norm(y-M(x,1))^2 ;
        f_previous = f_current;
        x = x + (1/mu)*M(y - M(x,1),2);
        x_tensor = reshape(x,[n1,n2,n3]) ;       
        [X,~,~] = prox_wtnn(x_tensor,hfun_sg,gamma,lambda,1/mu) ;        
        x = X(:) ;        
        f_current = norm(y-M(x,1)) + lambda*norm(x,1) ;
        if norm(f_current-f_previous)/norm(f_current + f_previous)<tol
            break;
        end
    end
    if norm(y-M(x,1))<err
        break ;
    end
    lambda = decfac*lambda;
end



