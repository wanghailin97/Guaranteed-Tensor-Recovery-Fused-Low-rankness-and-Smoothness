
% test on the tensor completion problem
currpath = cd ;
addpath(genpath(currpath)) ;


newDataFlag = 1 ;

if newDataFlag == 1    
    clc ;
    close all ;
    n = 100;
    n1 = n;
    n2 = n;
    n3 = n;
    r = 10; % tubal rank
    ML = randn(n1,r,n3)/n1; MR = randn(r,n2,n3)/n2;
    dr = (n1+n2-r)*r*n3;
    m = 3*dr;
    p = m/(n1*n2*n3);

   omega = find(rand(n1*n2*n3,1)<p);
end


IDX = omega ;
sizeX = [n1,n2,n3] ;
M = opRestriction(prod(sizeX), IDX);
X = tprod(ML,MR);
x = X(:);
y = M(x,1);

% IRNN with different penalties
% fun = 'lp' ;        gamma = 0.5 ;
 fun = 'scad' ;      gamma = 100 ;
% fun = 'logarithm' ; gamma = 10 ;
% fun = 'mcp' ; gamma = 10 ;
% fun = 'cappedl1' ; gamma = 1000 ;
% fun = 'etp' ;  gamma = 0.1 ;
% fun = 'geman' ;  gamma = 10 ;
% fun = 'laplace' ; gamma = 10 ;

tic
XRec = IRTNN(fun,y,M,n1,n2,n3,gamma) ;
toc

err = norm(X(:)-XRec(:),'fro')/norm(X(:),'fro')

tubalrank(XRec)





