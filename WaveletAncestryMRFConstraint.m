classdef WaveletAncestryMRFConstraint < handle
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
    end
    
    methods
        function constraint = WaveletAncestryMRFConstraint(texture,varargin)
            constraint.texture = texture;
            constraint.blocksize = 2;
            constraint.x_dataratio = .8;
            constraint.y_dataratio = .35;
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
            
            if ~isequal(constraint.x_sz,size(t.x))
                constraint.x_sz = size(t.x);
                constraint.update_exemplar_sample(); % exemplar has changed sizes: re sample
            end
            
            update_synthesis_sample(constraint); %resample patches from synthesis
            L = constraint.num_levels;
            x = t.x;
            y = t.y;
            x_swt0 = Transform.haar_udwt_2D(x,L);
            y_swt0 = Transform.haar_udwt_2D(y,L);
            C = size(t.y,3);
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
                X = Transform.im2col_sample(x_swt, constraint.blocksize, constraint.x_sample_ids);
                Y = Transform.im2col_sample(y_swt, constraint.blocksize, constraint.y_sample_ids);
                Y = Transport.nn_search(Y,X);
                y_swt = Transform.col2im_sample(Y,y_swt,constraint.blocksize,constraint.y_sample_ids);
                y_swt0(:,:,band_ids(1:3*C)) = y_swt(:,:,1:3*C);
            end
            t.y = Transform.haar_iudwt_2D(size(t.y),y_swt0);
        end
        
        function change_scale(constraint)
            update_exemplar_sample(constraint); 
        end
        
        function update_synthesis_sample(constraint)
            t = constraint.texture;
            constraint.x_sz = size(t.x);
            constraint.y_sample_ids = generate_sample_ids(size(t.y),constraint.blocksize,constraint.y_dataratio);
        end
        
        function update_exemplar_sample(constraint)
            t = constraint.texture;
            constraint.x_sample_ids = generate_sample_ids(size(t.x),constraint.blocksize,constraint.x_dataratio);
        end
    end
    
end

