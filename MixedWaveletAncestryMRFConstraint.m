classdef MixedWaveletAncestryMRFConstraint < handle
    properties (SetAccess = public)
        texture;
        blocksize;
        use_gpu = 0;
        x_sample_ids;
        y_sample_ids;
        x_dataratio;
        y_dataratio;
        x_sz; %x must be resampled if it changes size, e.g. a gaussian pyr
        num_levels
        wavelet_type %can only be dwt (haar), swt (haar), udtcwt (designed filter)
        
        band_weights; 
    end
    
    methods
        function constraint = MixedWaveletAncestryMRFConstraint(texture,varargin)
            constraint.texture = texture;
            constraint.blocksize = 2;
            constraint.x_dataratio = .3;
            constraint.y_dataratio = .3;
            
            if numel(varargin)
                params = varargin{1};
                if (isfield(params,'wavelet_type'));constraint.wavelet_type = params.wavelet_type;end;
                if (isfield(params,'num_levels'));constraint.num_levels = params.num_levels;end;
                if (isfield(params,'blocksize'));constraint.blocksize = params.blocksize;end;
                if (isfield(params,'x_dataratio'));constraint.x_dataratio = params.x_dataratio;end;
                if (isfield(params,'y_dataratio'));constraint.y_dataratio = params.y_dataratio;end;
            end
            constraint.x_sz =  0;
        end
        
        
        
        function enforce_constraint(constraint)
            t = constraint.texture;
            
            if ~isequal(constraint.x_sz,size(t.x1))
                constraint.x_sz = size(t.x1);
                constraint.update_exemplar_sample(); % exemplar has changed sizes: re sample
            end
            
            update_synthesis_sample(constraint); %resample patches from synthesis
            L = constraint.num_levels;
            x1 = t.x1;
            x2 = t.x2;
            y = t.y;
            x1_swt0 = Transform.haar_udwt_2D(x1,L);
            x2_swt0 = Transform.haar_udwt_2D(x2,L);
            y_swt0 = Transform.haar_udwt_2D(y,L);
            C = size(t.y,3);
            band_ids = [];
            for c = 1:C
               band_ids = [band_ids,c*(3*L+1)];
            end
            
            band_w  = constraint.texture.band_weights;
            for l = L:-1:1
                band_id = 3*(l-1)+1;
                prev_ids = band_ids;
                band_ids = [];
                ids = band_id:band_id+2;
                for c = 1:C
                    band_ids = [band_ids,ids+(c-1)*(3*L+1)];
                end
                band_ids = [band_ids,prev_ids];
                x1_swt = x1_swt0(:,:,band_ids);
                x2_swt = x2_swt0(:,:,band_ids);
                
                y_swt = y_swt0(:,:,band_ids);
                
                X = Transform.im2col_sample(x1_swt, constraint.blocksize, constraint.x_sample_ids);
                Y = Transform.im2col_sample(y_swt, constraint.blocksize, constraint.y_sample_ids);
                Y = Transport.nn_search(Y,X);
                y_swt1 = Transform.col2im_sample(Y,y_swt,constraint.blocksize,constraint.y_sample_ids);
                
                X = Transform.im2col_sample(x2_swt, constraint.blocksize, constraint.x_sample_ids);
                Y = Transform.im2col_sample(y_swt, constraint.blocksize, constraint.y_sample_ids);
                Y = Transport.nn_search(Y,X);
                y_swt2 = Transform.col2im_sample(Y,y_swt,constraint.blocksize,constraint.y_sample_ids);
                
                %blend according to band_weights
                for b = 1:3
                    y_swt0(:,:,band_ids(b:3:3*C)) = band_w(b)*y_swt1(:,:,b:3:3*C) + (1-band_w(b))*y_swt2(:,:,b:3:3*C);
                end
                
            end
            t.y = Transform.haar_iudwt_2D(size(t.y),y_swt0);
        end
        
        function x_swt = blend_swt(x1_swt,x2_swt,band_weights)

            num_levels = (size(x1_swt,3)-3)/9;

            band_weights = [repmat(band_weights(1:3),1,num_levels),band_weights(end)];%add lowpass assignment

            band_weights = repmat(band_weights,1,3);

            x_swt = zeros(size(x1_swt));

            for k = 1:size(x1_swt,3)
                x_swt(:,:,k) = linear_blend_images(x1_swt(:,:,k),x2_swt(:,:,k),band_weights(k));
            %     x_swt(:,:,k) = linear_blend_images_maxabs(x1_swt(:,:,k),x2_swt(:,:,k),band_weights(k));
            end
        end
        
        function change_scale(constraint)
            update_exemplar_sample(constraint); 
        end
        
        function update_synthesis_sample(constraint)
            t = constraint.texture;
            constraint.x_sz = size(t.x1);
            constraint.y_sample_ids = generate_sample_ids(size(t.y),constraint.y_dataratio);
        end
        
        function update_exemplar_sample(constraint)
            t = constraint.texture;
            constraint.x_sample_ids = generate_sample_ids(size(t.x1),constraint.x_dataratio);
        end
    end
    
end

