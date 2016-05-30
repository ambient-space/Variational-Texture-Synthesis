function [pos_map,copy_map] = get_tiling(x,y,blocksize)


[M,N,~] = size(x);
M = M-blocksize+1;N=N-blocksize+1;
[J,I] = meshgrid(1:N,1:M);
J = mat2gray(J); I = mat2gray(I);
pos_map = cat(3,J,I.*J,I);

ids_x = generate_sample_ids(size(x),blocksize,1);
ids_y = generate_sample_ids(size(y),blocksize,1);

[~,NN_ids] = NN_search_gpu(x,y,ids_x,ids_y,blocksize);

[M2,N2,~] = size(y);
M2 = M2-blocksize+1; N2 = N2 - blocksize+1;

copy_map = zeros(M2,N2,3);
for c = 1:3
    tmp1 = pos_map(:,:,c);
    copy_map(:,:,c) = reshape(tmp1(NN_ids),[M2,N2]);
end

r = imfilter(mean(copy_map,3),fspecial('laplacian'));
r = r>.015;

figure,
subplot(321);imshow(x);
subplot(322);imshow(y);
subplot(323);imshow(pos_map);
subplot(324);imshow(copy_map);
p = sum(r(:))/numel(r);
subplot(325);imshow(r); title(num2str(p));