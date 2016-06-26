classdef Analyze
    properties (Constant)
        
        
        
    end
    
    methods (Static)
        function synth_info = get_tiling(x,y,blocksize,show)
            
            [M,N,~] = size(x);
            pos_map = Analyze.get_pos_map(M,N);
            ids_x = generate_sample_ids(size(x),1);
            ids_y = generate_sample_ids(size(y),1);
            
            
            X = Transform.im2col_sample(x, blocksize, ids_x);
            Y = Transform.im2col_sample(y, blocksize, ids_y);
            
            [~,NN_ids] = Transport.nn_search(Y,X);
            [M2,N2,~] = size(y);
            
            tiling_image = zeros(M2,N2,3);
            for c = 1:3
                tmp1 = pos_map(:,:,c);
                tiling_image(:,:,c) = reshape(tmp1(NN_ids),[M2,N2]);
            end
            
            
            r = imfilter(mean(tiling_image,3),fspecial('laplacian'));
            r = r>.015;
            
            innovation_capacity = sum(r(:))/numel(r);
            if show
                figure,
                subplot(321);imshow(x);
                subplot(322);imshow(y);
                subplot(323);imshow(pos_map);
                subplot(324);imshow(tiling_image);
                
                subplot(325);imshow(r); title(num2str(innovation_capacity));
            end
            
            synth_info.innovation_capacity = innovation_capacity;
            synth_info.tiling_image = tiling_image;
            
        end
        
        function pos_map = get_pos_map(M,N)
            [J,I] = meshgrid(1:N,1:M);
            J = mat2gray(J); I = mat2gray(I);
            pos_map = cat(3,J,I.*J,I);
            
        end
    end
    
end

