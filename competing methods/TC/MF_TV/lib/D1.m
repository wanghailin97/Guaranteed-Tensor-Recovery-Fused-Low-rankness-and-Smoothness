function X=D1(X,n1,n2) %%% diff operator in vectical direction
for i=1:n2   %%ԭ�ȵ�
    temp=X(((i-1)*n1+1):(i*n1),:);
    X(((i-1)*n1+1):(i*n1),:)=[temp(2:end,:);temp(1,:)]-temp;
end

%%%�����������ǵȼ۵�,�ռ���Ĳ��
%frames = ; % ��Ƶ֡��
%for i=1:frames
%	D1 = Gradient(X(:,:,i),��ֱ����)��
%end

% X = [X(2:end,:);X(1,:)] - X; %%�Լ�д��