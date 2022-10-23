function A = inverselineartransform(A,L)

% Compute the result of  inverse linear  transform  L    on tensor A 
% Written by  Wenjin Qin  (qinwenjin2021@163.com)


p = length(size(A));
n = zeros(1,p);
for i = 1:p
    n(i) = size(A,i);
end


if strcmp(L,'fft')
     for i = p:-1:3
        A = ifft(A,[],i);
     end
elseif strcmp(L,'dct')
     for i = p:-1:3
        A = idct(A,[],i);
     end    
elseif iscell(L)
    for i = p:-1:3
        Mat_L=L{i-2};
        A = tmprod(A,inv(Mat_L),i);
    end
end