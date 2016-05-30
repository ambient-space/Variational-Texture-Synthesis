function crop_im = crop(source)

im = imread(source);


figure;

while 1
    imshow(source)
    disp('select center');
    [x,y]=ginput;
    M = input('enter image height: ');
    N = input('enter image width: ');
    B1 = floor(y-M/2:y+M/2-1);B2 = floor(x-N/2:x+N/2-1);
    if nnz(B1<1) || nnz(B2<1) || nnz(B1>size(im,1)) || nnz(B2>size(im,2))
        disp('invalid dimensions');
        continue;
    end
    im = im(B1,B2,:);
    imshow(im);
    if input('keep image? ')
        break;
    end
end

crop_im = im;