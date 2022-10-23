function A = lineartransform(A,L)

% Compute the result of  invertible  linear  transform  L    on tensor A 

% Written by  Wenjin Qin  (qinwenjin2021@163.com)




p = length(size(A));
n = zeros(1,p);
for i = 1:p
    n(i) = size(A,i);
end



if strcmp(L,'fft')
    for i = 3:p
        A = fft(A,[],i);
    end
elseif strcmp(L,'dct')
    for i = 3:p
        A = dct(A,[],i);
    end
elseif iscell(L)
    for i=1:(p-2)
        [l1,l2] = size(L{i});
        if l1 ~= l2 || l1 ~= n(i+2)
            error('Inner tensor dimensions must agree.');
        else
             A = tmprod(A,L{i},i+2);
        end
    end
   
end