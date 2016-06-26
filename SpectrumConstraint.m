classdef SpectrumConstraint < handle
    properties (SetAccess = public)
        texture
    end
    
    methods
        function constraint = SpectrumConstraint(texture)
            constraint.texture = texture;
        end
    end
    
    methods
        
        function enforce_constraint(constraint)
            x = constraint.texture.x;
            y = constraint.texture.y;
            
            constraint.texture.y = Spectrum.transfer(y,x);
            
        end
        
        
    end
end

