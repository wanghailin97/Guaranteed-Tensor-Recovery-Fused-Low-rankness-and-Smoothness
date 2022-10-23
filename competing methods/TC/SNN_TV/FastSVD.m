
%The function compute the left K singular vectors of matrix X,
%by using the probability method developed by Halko et al.

function [U,S,V] = FastSVD(X,K)
    %generate a random matrix with gaussian
     eta=2;
     [row, col]=size(X);
     Omega=randn(col,K);
     Y=(X*X')^eta*X*Omega;
     [Q,~]=qr(Y,0);
     B=Q'*X;
     [U,S,V]=svd(B,'econ');
     U=Q*U;
end

