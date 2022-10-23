%by Xutao Li
%The code is for the LRTC_TV_I method in 2017 AAAI conference paper "Low Rank Tensor Completeion with Total Variation for Visual Data Inpainting"



% LRTC_TV_I is used for recover a tensor from some observed
% entries, with Tensor unfolding formulation

% lambda is a scalar used as parameter for TV terms
% alpha  is a vector used as parameter for tensor trace norm
% beta is a binary vector used as parameter to indicate constraints or not
% index is used to record the set of indices of observations
% value records the values of observations
% N is the order of the tensor
% tsize is the size of the tensor
function [ tensor_Z ] = LRTC_TV_I(index, value, lambda, alpha, beta, tsize, N )
    
    M=cell(N,1);
    Q=cell(N,1);
    F=cell(N,1);
    Z=cell(N,1);
    R=cell(N,1);
    Lambda=cell(N,1);
    Gamma=cell(N,1);
    Phi=cell(N,1);
    
%    epsilon=1.0e-6;
     epsilon=1.0e-8;
    for i=1:N
        l=1;
        for j=[1:i-1,i+1:N]
            l=l*tsize(j);
        end
        M{i}=zeros(tsize(i),l);
        if(beta(i)==1)
            Q{i}=rand(tsize(i)-1, l);
            R{i}=rand(tsize(i),l);
        else
            Q{i}=[];
            R{i}=[];
        end
       
        F{i}=zeros(tsize(i)-1,tsize(i));
        for j=1:tsize(i)-1
            F{i}(j,j)=1;
            F{i}(j,j+1)=-1;
        end

    end
    for i=1:N
       if(beta(i)==1)
           Lambda{i}=sign(Q{i})/max([norm(Q{i}), norm(Q{i},Inf), epsilon]);
           Phi{i}=sign(R{i})/max([norm(R{i}), norm(R{i},Inf), epsilon]);
       else
           Lambda{i}=[];
           Phi{i}=[];
       end    
       Gamma{i}=sign(M{i})/max([norm(M{i}), norm(M{i},Inf), epsilon]);
       
    end
    
    tensor_Z=tenzeros(tsize);
    tensor_Z(index)=value;
    
    iteration=1;
%    myInitial_v=1.0e-3;
    myInitial_v=1.0e-4;
    rho=myInitial_v;
    mu=myInitial_v;
    gamma=myInitial_v;
%   factor=1.05;
    factor = 1.1; 
while(true)
        tensor_Z_pre=tensor_Z;
        for n=1:N
            Z{n}=double(tenmat(tensor_Z,n));
        end

        %update Q
        for n=1:N
            if(beta(n)==1)
                Q{n}=myshrinkage(F{n}*R{n} - 1/rho *Lambda{n}, lambda/rho);
            end
        end
        %update M
        for n=1:N
            tmpMatrix=Z{n}  -  1/mu*Gamma{n};
            M{n}=SVT(tmpMatrix, alpha(n)/mu,0);
        end
        %update R
        for n=1:N
            if(beta(n)==1)
                R{n}= (rho*F{n}'*F{n} + gamma * eye((tsize(n))))\(F{n}'* Lambda{n}+ rho* F{n}'* Q{n} + gamma * Z{n} - Phi{n});
            end
        end
        
        
        %update Z
        tmp=0;
        for n=1:N
            if(beta(n)==1)
               currZ= Gamma{n} + Phi{n} +  mu * M{n} + gamma * R{n};
            else
               currZ= Gamma{n} +  mu * M{n};
            end
            myZ=tenmat(tensor_Z,n);
            myZ=tenmat(currZ,myZ.rdims, myZ.cdims,myZ.tsize);
            tmp=tmp + tensor(myZ);
        end
        NN=numel(find(beta==1));
        
        tensor_Z=tmp/(N*mu + NN*gamma);
        tensor_Z(index)=value;
        
        for n=1:N
            if(beta(n)==1)
                Lambda{n}=Lambda{n} + rho * (Q{n}- F{n}*R{n});
                Phi{n}= Phi{n}+ gamma*(R{n} - double(tenmat(tensor_Z,n)));
            end
            Gamma{n}=Gamma{n}+ mu*(M{n} - double(tenmat(tensor_Z,n)));
        end
        
        diff=norm(tensor_Z-tensor_Z_pre)/norm(tensor_Z);
        
        larange_cond_1=0;
        larange_cond_2=0;
        larange_cond_3=0;
        for n=1:N
            if(beta(n)==1)
                larange_cond_1= larange_cond_1 +  sum(sum(Lambda{n}.* (Q{n} - F{n}*R{n})));
                larange_cond_3= larange_cond_3 +  sum(sum(Phi{n}.*(R{n} - double(tenmat(tensor_Z,n)))));
            end
            larange_cond_2= larange_cond_2 +  sum(sum(Gamma{n}.* (M{n} - double(tenmat(tensor_Z,n)))));
        end

        
        rho=rho*factor;
        mu=mu*factor;
        gamma=gamma*factor;
        
        fprintf('iter=%d,diff=%f\n',iteration,diff);
        fprintf('iter=%d,lagrange_cond=%f,lagrange_cond_2=%f,lagrange_cond_3=%f\n',iteration,larange_cond_1,larange_cond_2,larange_cond_3);
        
        if(diff<epsilon||iteration>500)
            break;
        end
        
        %figure(4);imshow(double(tensor_Z));
        iteration = iteration + 1;
        
    end
    tensor_Z=double(tensor_Z);
end

