function X=D2(X,n1,n2) %%% diff operator in hor direction
X=[X((n1+1):n1*n2,:);X(1:n1,:)]-X;   %%ԭ�ȵ�

% X = [X(:,2:end),X(:,1)] - X; %%�Լ�д��