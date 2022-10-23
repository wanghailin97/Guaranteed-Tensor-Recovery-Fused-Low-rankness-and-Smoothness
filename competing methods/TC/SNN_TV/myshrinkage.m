function [A]=myshrinkage(A,a)
    [row, col]=size(A);
    for i=1:row
        for j=1:col
              if(A(i,j)==0)
                  continue;
              end
              if(abs(A(i,j))<a)
                  A(i,j)=0;
              else
                  A(i,j)=A(i,j)*(1-a/abs(A(i,j)));
              end
        end
    end
end