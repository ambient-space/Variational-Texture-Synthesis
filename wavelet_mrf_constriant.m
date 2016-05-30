function synthesis = wavelet_mrf_constriant(synthesis, exemplar, blocksize, s_sample_ratio, e_sample_ids,L)
%INPUTS:
%synthesis = current RGB synthesis image




s_sample_ids = generate_sample_ids(size(synthesis),blocksize,s_sample_ratio);
x_swt0 = Transform.haar_udwt_2D(exemplar,L);
y_swt0 = Transform.haar_udwt_2D(synthesis,L);
C = size(synthesis,3);
band_ids = [];
for c = 1:C
    band_ids = [band_ids,c*(3*L+1)];
end

for l = L:-1:1
    band_id = 3*(l-1)+1;
    prev_ids = band_ids;
    band_ids = [];
    ids = band_id:band_id+2;
    for c = 1:C
        band_ids = [band_ids,ids+(c-1)*(3*L+1)];
    end
    band_ids = [band_ids,prev_ids];
    x_swt = x_swt0(:,:,band_ids);
    y_swt = y_swt0(:,:,band_ids);
    X = Transform.im2col_sample(x_swt, blocksize, e_sample_ids);
    Y = Transform.im2col_sample(y_swt, blocksize, s_sample_ids);
    Y = Transport.nn_search(Y,X);
    y_swt = Transform.col2im_sample(Y,y_swt,blocksize,s_sample_ids);
    y_swt0(:,:,band_ids(1:3*C)) = y_swt(:,:,1:3*C);
end
synthesis = Transform.haar_iudwt_2D(size(synthesis),y_swt0);
