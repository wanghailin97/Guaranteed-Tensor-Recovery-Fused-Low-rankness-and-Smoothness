function [result] = rse(X_ground,Y)
    X=double(X_ground);
    Y=double(Y);
    a=(X(:)-Y(:)).^2;
    b=(X(:)).^2;
    a=sum(a(:));
    b=sum(b(:));
    result=sqrt(a)/sqrt(b);
end

