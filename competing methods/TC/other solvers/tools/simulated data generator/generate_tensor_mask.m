function coef_tensor=generate_tensor_mask(row_1,row_2,rk)
mask = zeros(row_1,row_2);
coef_tensor = zeros(row_1,row_2,rk);
class = rk;
support = normrnd(0,1,class,rk);
%support = rand(class,rk);
init_center = randperm(row_1*row_2,class);
center_axis = zeros(class,2);
for j = 1:class
    row_id = mod(init_center(j),row_1);
    col_id = (init_center(j)-row_id)/row_1+1;
    center_axis(j,1) = row_id;
    center_axis(j,2) = col_id;
end
for j = 1:row_2
    for i = 1:row_1
        a = repmat([j,i],class,1);
        dist = sum((a-center_axis).*(a-center_axis),2);
        [c,c_id]=min(dist);
        mask(j,i)=c_id;
        coef_tensor(j,i,:)=support(c_id,:);
    end
end
%coef_matrix = reshape(coef_tensor,[row_1*row_2,rk]);


            