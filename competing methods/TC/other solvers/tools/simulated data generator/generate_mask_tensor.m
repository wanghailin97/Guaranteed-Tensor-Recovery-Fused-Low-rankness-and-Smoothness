function [data,endmember] = generate_mask_tensor(rank_number,step)
%% generate endmember matrix
endmember = normrnd(0,1,rank_number,200);
for i=1:rank_number
    endmember(i,:) = smooth(endmember(i,:));
end

[mask,total_coef_matrix]=generate_mask(10*step,10*step,rank_number);
data = total_coef_matrix*endmember;
for i=1:200
    data(:,i)= (data(:,i)-min(data(:,i)))/(max(data(:,i))-min(data(:,i))+0.0001);
end

        


