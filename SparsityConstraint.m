classdef SparsityConstraint < handle
    properties (SetAccess = public)
        texture;
        blocksize;
        
        y_sample_ids;
        x_dataratio;
        y_dataratio;
        x_sz; %x must be resampled if it changes size, e.g. a gaussian pyr
        
        dictionaries;
        sparsity;
        dictsize;
    end
    
    methods
        function constraint = SparsityConstraint(texture,varargin)
            constraint.texture = texture;
            %Default properties
            constraint.blocksize = 10;
            constraint.x_dataratio = .4;
            constraint.y_dataratio = .5;
            constraint.sparsity = 3;
            constraint.dictsize = 2^8;
            
            %Read in specified properties
            if numel(varargin)
                params = varargin{1};
                read_params(constraint,params);
            end
            
            constraint.x_sz =  0;
            constraint.dictionaries =  {};
        end
        
        
        
    end
    
    methods
        
        function enforce_constraint(constraint)
            t = constraint.texture;
            
            if ~isequal(constraint.x_sz,size(t.x))
                constraint.x_sz = size(t.x);
                constraint.update_dictionary();
            end
            
            update_synthesis_sample(constraint); %resample patches from synthesis
            
            y = t.y;
            D = constraint.dictionaries{numel(constraint.dictionaries)}; % get most recent dict
            
            Y = Transform.im2col_sample(y,constraint.blocksize,constraint.y_sample_ids);

            G = Transport.sparsecode(Y,D,constraint.sparsity);

            Y = D*G;
            
            t.y =  Transform.col2im_sample(Y,y,constraint.blocksize,constraint.y_sample_ids);
            
        end
        
        function change_scale(constraint)
            update_exemplar_sample(constraint); 
        end
        
        function update_synthesis_sample(constraint)
            t = constraint.texture;
            constraint.x_sz = size(t.x);
            constraint.y_sample_ids = generate_sample_ids(size(t.y),constraint.y_dataratio);
        end
        
        function update_dictionary(constraint)
            t = constraint.texture;
            
            x_sample_ids = generate_sample_ids(size(t.x),constraint.x_dataratio);
            
            X = Transform.im2col_sample(t.x,constraint.blocksize,x_sample_ids);
            
            [D,~] = ksvd(X,constraint.dictsize,constraint.sparsity);
            constraint.dictionaries{numel(constraint.dictionaries)+1} = D;
           
        end
        
    end
end

