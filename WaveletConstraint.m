classdef WaveletConstraint < handle
    properties (SetAccess = public)
        texture
        num_levels
        wavelet_type %can only be dwt (haar), swt (haar), udtcwt (designed filter)
        
    end
    
    methods
        function constraint = WaveletConstraint(texture,varargin)
            constraint.texture = texture;
            constraint.num_levels = 4;
            constraint.wavelet_type = 'haar_udwt';
            if numel(varargin)
                params = varargin{1};
                if (isfield(params,'wavelet_type'));constraint.wavelet_type = params.wavelet_type;end;
                if (isfield(params,'num_levels'));constraint.num_levels = params.num_levels;end;
                if (isfield(params,'scales'));constraint.num_levels = constraint.num_levels;end;% - numel(params.scales);end;
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
        
        function change_scale(constraint)
%             constraint.num_levels = constraint.num_levels + 1;
        end
        
    end
end

