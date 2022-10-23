function Ten_mask=generate_tensor_mask_4D(n1,n2,n3,n4, k)
Ten_mask = zeros(n1, n2, n3, n4);
Ten_mask_value = normrnd(0,1,k,1);
center = randperm(n1*n2*n3*n4,k);
center_axis = zeros(k,4);
for i = 1:k
    [center_axis(i,1), center_axis(i,2),center_axis(i,3),center_axis(i,4)]= ind2sub([n1,n2,n3,n4],center(i));
end
for axis_1 = 1:n1
    for axis_2 = 1:n2
        for axis_3 = 1:n3
            for axis_4 = 1:n4
                a = repmat([axis_1, axis_2, axis_3, axis_4], k, 1);
                dist = sum((a-center_axis).*(a-center_axis),2);
                [~, close_id] = min(dist);
                Ten_mask(axis_1, axis_2, axis_3, axis_4) = Ten_mask_value(close_id);
            end
        end
    end
end



            