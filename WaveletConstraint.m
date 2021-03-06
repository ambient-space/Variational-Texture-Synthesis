classdef WaveletConstraint < handle
    properties (SetAccess = public)
        texture
        num_levels
        wavelet_type %can only be a haar udwt for now
        
    end
    
    methods
        function constraint = WaveletConstraint(texture,varargin)
            constraint.texture = texture;
            constraint.num_levels = 4;
            constraint.wavelet_type = 'haar_udwt';
            if numel(varargin)
                params = varargin{1};
                read_params(constraint,params);
            end
        end
    end
    
    methods
        
        function enforce_constraint(constraint)
            x = constraint.texture.x;
            y = constraint.texture.y;
            L = constraint.num_levels;
            
            [x_pca,base] = Transform.pca(x);
            [y_pca,~] = Transform.pca(y,base);
            x_udwt = Transform.haar_udwt_2D(x_pca,L);
            y_udwt = Transform.haar_udwt_2D(y_pca,L);
            for b = 1:size(x_udwt,3)
                y_udwt(:,:,b) = Transport.optimal_transport_1D(y_udwt(:,:,b),x_udwt(:,:,b));
            end
            y_pca = Transform.haar_iudwt_2D(size(y),y_udwt);
            y = Transform.ipca(y_pca,base);
            
            constraint.texture.y = y;
        end
        
    end
end

