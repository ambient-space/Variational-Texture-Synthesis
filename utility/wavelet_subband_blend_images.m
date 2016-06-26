function x = wavelet_subband_blend_images(x1,x2,band_weights,num_levels)
%band_weights = weights over subbands transform hh,lh,hl,ll

x1_swt = Transform.haar_udwt_2D(x1,num_levels);
x2_swt = Transform.haar_udwt_2D(x2,num_levels);

band_weights = [repmat(band_weights(1:3),1,num_levels),band_weights(end)];%add lowpass assignment

band_weights = repmat(band_weights,1,3);

x_swt = zeros(size(x1_swt));

for k = 1:size(x1_swt,3)
    x_swt(:,:,k) = band_weights(k)*x1_swt(:,:,k) +(1-band_weights(k))* x2_swt(:,:,k);
end


x = Transform.haar_iudwt_2D(size(x1),x_swt);