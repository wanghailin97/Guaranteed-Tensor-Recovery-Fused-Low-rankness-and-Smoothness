function Mat_mask=generate_matrix_mask_new(n1,n2,k)
Mat_mask = zeros(n1, n2);
Mat_mask_value = normrnd(0,1,k,1);
center = randperm(n1*n2,k);
center_axis = zeros(k,2);
for i = 1:k
    [center_axis(i,1), center_axis(i,2)]= ind2sub([n1, n2],center(i));
end
for axis_1 = 1:n1
    for axis_2 = 1:n2
            a = repmat([axis_1, axis_2], k, 1);
            dist = sum((a-center_axis).*(a-center_axis),2);
            [~, close_id] = min(dist);
            Mat_mask(axis_1, axis_2) = Mat_mask_value(close_id);
    end
end



            