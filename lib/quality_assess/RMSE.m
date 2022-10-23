function rmse = RMSE(hsi1, hsi2)
rmse = sqrt(sum((hsi1(:)-hsi2(:)).^2)/numel(hsi1));
end


