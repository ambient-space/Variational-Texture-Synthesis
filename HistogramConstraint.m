classdef HistogramConstraint < handle
    properties (SetAccess = public)
        texture
    end
    
    methods
        function constraint = HistogramConstraint(texture)
            constraint.texture = texture;
        end
    end
    
    methods
        
        function enforce_constraint(constraint)
            x = constraint.texture.x;
            y = constraint.texture.y;
            [x_pca,base] = Transform.pca(x);
            [y_pca,~] = Transform.pca(y,base);
            for c  = 1:3
                y_pca(:,:,c) = Transport.optimal_transport_1D(y_pca(:,:,c),x_pca(:,:,c));
            end
            y = Transform.ipca(y_pca,base);
            constraint.texture.y = y;
        end
        
    end
end

