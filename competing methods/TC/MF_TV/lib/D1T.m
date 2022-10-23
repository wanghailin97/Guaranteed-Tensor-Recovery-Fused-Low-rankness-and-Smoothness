function X=D1T(X,n1,n2) %%%transpose of diff operator in vectical direction
for i=1:n2  %%原先的
    temp=X(((i-1)*n1+1):(i*n1),:);
    X(((i-1)*n1+1):(i*n1),:)=[temp(end,:);temp(1:end-1,:)]-temp;
end

% X = [X(end,:);X(1:end-1,:)] - X; %%自己写的