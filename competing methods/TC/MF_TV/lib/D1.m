function X=D1(X,n1,n2) %%% diff operator in vectical direction
for i=1:n2   %%原先的
    temp=X(((i-1)*n1+1):(i*n1),:);
    X(((i-1)*n1+1):(i*n1),:)=[temp(2:end,:);temp(1,:)]-temp;
end

%%%以上与下面是等价的,空间域的差分
%frames = ; % 视频帧数
%for i=1:frames
%	D1 = Gradient(X(:,:,i),竖直方向)；
%end

% X = [X(2:end,:);X(1,:)] - X; %%自己写的