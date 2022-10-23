function [X,iter,Res,Rel]=TV_X3(Y,A,X_true,mu,beta,rho,miter,initer,n1,n2)
%%
% sloving X3 subproblem.
% using SUnADM with isotropic TV

% created by Teng-Yu Ji 
% 1/13/2018

%%
%%
m=size(A,1);
N=size(Y,1);
%%%eigenvalues of AAT
[U1,S1,V1]=svd(A*A');
Sig1=diag(S1);
%%%eigenvalues of  BTB=D1TD1+D2TD2
d1=zeros(n1,n2);
d1(end,1)=1;d1(1,1)=-1;
d2=zeros(n1,n2);
d2(1,end)=1;d2(1,1)=-1;
d1eig=fft2(d1);
d2eig=fft2(d2);
Sig2=beta*((abs(d1eig)).^2+(abs(d2eig)).^2);
%%%eigenvalues
Sig=repmat(Sig1',N,1)+repmat(Sig2(:),1,m)+rho;
Sig=1./Sig;
%%% initialization
Lam=zeros(2*N,m);
W1=zeros(N,m);
W2=zeros(N,m);
X_p=zeros(size(X_true));
X_p(1,1)=1;
% X_p=X_true;
res=1;
Res=[];
SRE=[];
Rel=[];
tol=10^(-3);
iter=0;
while res>tol && iter<miter
    for j=1:initer
    %%%X subproblem by solving Sylvester matrix euqation
    M=Y*A'+beta*(D1T(W1,n1,n2)+D2T(W2,n1,n2))-(D1T(Lam(1:N,:),n1,n2)+D2T(Lam((N+1):2*N,:),n1,n2))+rho*X_true;
    temp=Sig.*(calF(M)*U1);
    X=real(calinvF(temp))*U1';
    %%W subproblem
    Z1 = D1(X,n1,n2) + Lam(1:N,:)/beta;
    Z2 = D2(X,n1,n2) + Lam((N+1):2*N,:)/beta;
    V1=Z1.^2;
    V2=Z2.^2;
    sV=sqrt(V1+V2);
    sV(sV==0) = 1;
    sV= max(sV - mu/beta, 0)./sV;
    W1 = Z1.*sV;
    W2 = Z2.*sV;
    %caculated relative errors
    res=norm(X-X_p,'fro')/norm(X_p,'fro');
    Res=[Res,res];
    X_p=X;
    end
    %%% updating Lam  
    Lam=Lam+beta*([D1(X,n1,n2);D2(X,n1,n2)]-[W1;W2]);
    Rel=[Rel,norm(X-X_true,'fro')/norm(X_true,'fro')];
    iter=iter+1;
end
% X(X<0)=0;
    function X=calF(X) %%%calculate FX
        [N,m]=size(X);
        for ii=1:m
            temp1=fft2(reshape(X(:,ii),n1,n2));
            X(:,ii)=reshape(temp1,N,1);
        end
    end
    function X=calinvF(X) %%%calculate F^{*}X
        [N,m]=size(X);
        for ii=1:m
            temp2=ifft2(reshape(X(:,ii),n1,n2));
            X(:,ii)=reshape(temp2,N,1);
        end
    end
end