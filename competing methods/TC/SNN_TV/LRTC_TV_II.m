%by Xutao Li
%The code is for the LRTC_TV_II method in 2017 AAAI conference paper "Low Rank Tensor Completeion with Total Variation for Visual Data Inpainting"


% LRTC_TV_II is used for recover a tensor from some observed
% entries, with Tensor decomposition formulation

% lambda is a scalar used as parameter for TV terms
% alpha  is a vector used as parameter for tensor trace norm
% beta is a binary vector used as parameter to indicate constraints or not
% index is used to record the set of indices of observations
% value records the values of observations
% N is the order of the tensor
% tsize is the size of the tensor


function [tensor_Z  ] = LRTC_TV_II(index, value, lambda_1,lambda_2, alpha, beta, tsize, N )
    Q=cell(N,1);
    F=cell(N,1);
    R=cell(N,1);
    U=cell(N,1);
    V=cell(N,1);
    myV=cell(N,1);
    
    G=cell(N,1);
    Z=cell(N,1);
    W=cell(N,1);
    
    Lambda=cell(N,1);
    Gamma=cell(N,1);
    Phi=cell(N,1);

    
    epsilon=1.0e-6;
    for i=1:N
        l=1;
        for j=[1:i-1,i+1:N]
            l=l*tsize(j);
        end
        U{i}=rand(tsize(i),tsize(i));
        V{i}=rand(tsize(i),tsize(i));
       % V{i}=V{i}/norm(V{i});
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
       Gamma{i}=sign(U{i})/max([norm(U{i}), norm(U{i},Inf), epsilon]);
       
    end
    
    tensor_Z=tenzeros(tsize);
    tensor_Z(index)=value;
    tensor_G=tensor_Z;
    tensor_W=tenzeros(tsize);
    
    
    iteration=1;
    
    myInitial_v=1.0e-2;
    rho_1=myInitial_v;
    rho_2=myInitial_v;
    rho_3=myInitial_v;
    rho_4=myInitial_v;
    
    factor=1.1;
    
    while(true)
        tensor_Z_pre=tensor_Z;
        for n=1:N
            Z{n}=double(tenmat(tensor_Z,n));
            W{n}=double(tenmat(tensor_W,n));
        end
        
        %update Q
        for n=1:N
            if(beta(n)==1)
                Q{n}=myshrinkage(F{n}*R{n} - 1/rho_1 *Lambda{n}, lambda_1/rho_1);
            end
        end
        
        %update R
        for n=1:N
            if(beta(n)==1)
                R{n}= (rho_1*F{n}'*F{n} + rho_2 * eye((tsize(n))))\(F{n}'* Lambda{n}+ rho_1* F{n}'* Q{n} + rho_2 * Z{n} - Phi{n});
            end
        end
        
        %update U
        for n=1:N
            U{n}=SVT(V{n} + Gamma{n}/rho_3,alpha(n)/rho_3,0);
        end
        
        %update V
        for n=1:N
            tmp=ttm(tensor_G,V,-n);
            tmp=tenmat(tmp,n);
            tmp=double(tmp);
            V{n}=(- Gamma{n} + rho_3 * U{n} + W{n} * tmp' + rho_4 * Z{n}* tmp')/((rho_3*eye(tsize(n)) + rho_4*(tmp*tmp'))); 
        end
        
        %update Z
        tmp=0;
        for n=1:N
            if(beta(n)==1)
               currZ = Phi{n} + rho_2 * R{n};
               myZ=tenmat(tensor_Z,n);
               myZ=tenmat(currZ,myZ.rdims, myZ.cdims,myZ.tsize);
               tmp=tmp + tensor(myZ);
             end
        end
        tmp = tmp - tensor_W + rho_4* ttm(tensor_G,V,1:N);
        
        NN=numel(find(beta==1));
        tensor_Z = tmp/(NN*rho_2 + rho_4);
        tensor_Z(index)=value;
        
        %update G
%         n=3;
%         V_negative_n=1;
%         V_positive_n=V{n};
%         myDims=1;
%         for i=[1:n-1, n+1,N]
%             V_negative_n=kron(V_negative_n, V{i});
%             myDims=myDims*tsize(i);
%         end
%         tmp1=V_negative_n'*V_negative_n;
%         tmp2=V_positive_n'*V_positive_n;
%         tmp=kron(tmp1,tmp2);
%         tmp=tmp + lambda_2 * eye(size(tmp));
%         
%         
%         tmp1 = V_positive_n'* W{n}* V_positive_n;
%         tmp2 = rho_4 * V_positive_n'* Z{n} * V_positive_n;
%         
%         myG = tenmat(tensor_G, n);
%         myG_vector= tmp\(tmp1+tmp2);
%         myG_matrix=reshape(myG_vector,n,myDims);
%         myG=tenmat(myG_matrix,myG.rdims, myG.cdims,myG.tsize);
%         tensor_G=tensor(myG);

        for n=1:N
            myV{n}=V{n}';
        end
         myG = optimize_Z(myV,double(tensor_Z),double(tensor_W),rho_4,lambda_2);
         tensor_G = tensor(myG);
        
        
        %update multiplers
        for n=1:N
            if(beta(n)==1)
                Lambda{n}=Lambda{n} + rho_1 * (Q{n}- F{n}*R{n});
                Phi{n}= Phi{n}+ rho_2*(R{n} - double(tenmat(tensor_Z,n)));
            end
            Gamma{n}=Gamma{n}+ rho_3*(V{n}-U{n});
        end
        tensor_W = tensor_W + rho_4*(tensor_Z - ttm(tensor_G, V,1:N));
        diff=norm(tensor_Z-tensor_Z_pre)/norm(tensor_Z);
        
        larange_cond_1=0;
        larange_cond_2=0;
        larange_cond_3=0;
        larange_cond_4=0;
        for n=1:N
            if(beta(n)==1)
                larange_cond_1= larange_cond_1 +  sum(sum(Lambda{n}.* (Q{n} - F{n}*R{n})));
                larange_cond_2= larange_cond_2 +  sum(sum(Phi{n}.*(R{n} - double(tenmat(tensor_Z,n)))));
            end
            larange_cond_3= larange_cond_3 +  sum(sum(Gamma{n}.* (V{n} - U{n})));
        end
        tmp= tensor_W .* (tensor_Z - ttm(tensor_G, V,1:N));
        tmp=double(tmp);
        for i=1:tsize(1)
            for j=1:tsize(2)
                for k=1:tsize(3)
                    larange_cond_4 = larange_cond_4 + tmp(i,j,k);
                end
            end
        end
        
        rho_1=rho_1*factor;
        rho_2=rho_2*factor;
        rho_3=rho_3*factor;
        rho_4=rho_4*factor;
        
        fprintf('iter=%d,diff=%f\n',iteration,diff);
        fprintf('iter=%d,lagrange_cond=%f,lagrange_cond_2=%f,lagrange_cond_3=%f,,lagrange_cond_4=%f\n',iteration,larange_cond_1,larange_cond_2,larange_cond_3,larange_cond_4);
        
        if(iteration>300)
            break;
        end
        
        %double(tensor_Z)
        figure(9);imshow(double(tensor_Z));
        iteration = iteration + 1;
                
    end
    tensor_Z=double(tensor_Z);

end

