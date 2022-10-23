function X=D2T(X,n1,n2) %%%transpose of diff operator in horizatal direction
X=[X((n1*(n2-1)+1):n1*n2,:);X(1:(n1*(n2-1)),:)]-X;   %%原先的

% X = [X(:,end),X(:,1:end-1)] - X; %%自己写的