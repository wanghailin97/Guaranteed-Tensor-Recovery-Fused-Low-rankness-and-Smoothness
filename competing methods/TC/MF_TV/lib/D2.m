function X=D2(X,n1,n2) %%% diff operator in hor direction
X=[X((n1+1):n1*n2,:);X(1:n1,:)]-X;   %%原先的

% X = [X(:,2:end),X(:,1)] - X; %%自己写的