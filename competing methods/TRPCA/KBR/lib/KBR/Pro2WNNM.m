function [X, n, SigmaNew,wNorm] = Pro2WNNM(Z, tau)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% min: 1/2*||Z-X||^2 + ||X||_tr
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% [S, V, D, Sigma2] = MySVDtau(Z, tau);
% V = max(diag(V) - tau, 0);
% n = sum(V > 0);
% X = S(:, 1:n) * diag(V(1:n)) * D(:, 1:n)';

%% new
[m, n] = size(Z);
if m < n
    AAT = Z*Z';
    [S, Sigma, D] = svd(AAT);
    Sigma         = sqrt(diag(Sigma));    
    tol           = max(size(Z)) * eps(max(Sigma));
    [SigmaNew,n]  = ClosedWNNM(Sigma,tau,tol);
    wNorm         = sum(SigmaNew./(Sigma(1:n)+tol));
    SigmaNew      = SigmaNew ./ Sigma(1:n);
    X = S(:, 1:n) * diag(SigmaNew) * S(:, 1:n)' * Z;
    return;
end
if m > n
    [X, n, SigmaNew, wNorm] = Pro2WNNM(Z', tau);
    X = X';
    return;
end

