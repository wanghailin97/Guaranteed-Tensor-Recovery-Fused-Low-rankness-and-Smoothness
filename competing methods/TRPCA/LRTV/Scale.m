function [scaled]=Scale(data,lower,upper)
if(nargin <3)
    lower=-1;
    upper=1;
elseif(lower>upper)
    disp(['inputting error'])
end
[M,N]=size(data);
maxV=max(max(data));
minV=min(min(data));
scaled= (upper-lower)*(data-ones(M,N)*minV)./(maxV-minV);