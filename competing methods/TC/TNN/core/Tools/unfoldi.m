function I = unfoldi(i,j,L)

% Written by  Wenjin Qin  (qinwenjin2021@163.com)

I = ones(1,j);
    for t = j:-1:3   
        I(t) = ceil(i/L(t-1));   
        i = i-(I(t)-1)*L(t-1);
    end
end
