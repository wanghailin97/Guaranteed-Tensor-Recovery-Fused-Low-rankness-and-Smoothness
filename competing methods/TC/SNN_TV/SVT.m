
function [ myM ] = SVT(myM, mu,speedup)
    if(speedup==0)
        [U,S,V]=svd(myM,'econ');
    else
        [U,S,V]=FastSVD(myM,100);
    end
    [row,col]=size(S);
    bound=min([row col]);
    for i=1:bound
        if(S(i,i)>mu)
            S(i,i)=S(i,i)-mu;
        else
            S(i,i)=0;
        end
    end
    myM=U*S*V';
end

