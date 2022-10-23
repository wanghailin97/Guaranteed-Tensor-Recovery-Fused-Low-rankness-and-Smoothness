function [ Img ] = PseudoImage(hsi, slected_bands)
hsi_sz  =  size(hsi);
Img = zeros(hsi_sz(1), hsi_sz(2), 3);
Img(:,:,1) = hsi(:,:,slected_bands(1));
Img(:,:,2) = hsi(:,:,slected_bands(2));
Img(:,:,3) = hsi(:,:,slected_bands(3));
end




