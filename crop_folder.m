
addpath('./arbyreed raw textures/');
files = dir('./arbyreed raw textures/*.jpg');

formatted_folder = './arbyreed textures/';
mkdir(formatted_folder);
figure;
out_size = [480,640];
r = out_size(1)/out_size(2);
for file = files'
    
    im = imread(file.name);
    imshow(im)
    
    while 1
        imshow(im)
        disp('top left of bbox, btm right of bbox');
        
        [x1,y1] = ginput;
        [x2,y2] = ginput;
        
        
        H = y2-y1+1;
        W = x2-x1+1;
        
        if H/W > r
            
            %height is disproportional
            H = ceil( W*r);
            
        elseif H/W < r
            
            W = ceil( H/r);
            
        end
        
        B1 = max( floor( ( y1: min(y1+H-1,size(im,1)) ) ) , 1);
        B2 = max( floor( ( x1: min(x1+W-1,size(im,2)) ) ) , 1);
        
        tmp0 = im(B1,B2,:);
        tmp = imresize(im(B1,B2,:),out_size);
        imshow(tmp);
        if input('keep image? ')
            
            imwrite(tmp,[formatted_folder,file.name(1:end-4),'.png'],'png');
            
            break;
        end
    end
    
    crop_im = im;
    
    
end