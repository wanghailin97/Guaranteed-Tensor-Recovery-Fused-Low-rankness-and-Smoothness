function i = foldi(n,j,L)

% Written by  Wenjin Qin  (qinwenjin2021@163.com)

i = 1;
for t = j:-1:3   
   i = i + (n(t)-1)*L(t-1);
end
