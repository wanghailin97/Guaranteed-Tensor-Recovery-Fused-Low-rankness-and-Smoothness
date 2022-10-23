function N =nc(I,j,n)

% Written by  Wenjin Qin  (qinwenjin2021@163.com)

N = zeros(1,length(n));
if j >= 3
    for t = 3 : j-1
        if I(t) == 1 
            N(t) = 1;
        else
            N(t) = n(t)-I(t)+2;
        end
    end
end

N(j) = n(j)-I(j)+2;

if j<length(n)
for t = j+1 : length(n)
    N(t) = 1;
end
end 

end