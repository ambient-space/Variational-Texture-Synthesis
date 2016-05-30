classdef Texture < handle
    
    properties (SetAccess = public)
        x0
        y0
        iter
        scales
        in_scale
        out_scale
        type
        constraints = {};
        x
        y
    end
    
    
    
    methods
        function texture = Texture(x_in,varargin)
            texture.x0 = x_in;
            %default parameters
            texture.y0 =  rand(size(texture.x0));
            texture.iter = 8;
            texture.scales = [1,.5,.25];
            texture.in_scale = 1;
            texture.out_scale = 1;
            texture.type = 'alternating';
            
            if numel(varargin)
                params = varargin{1};
                if (isfield(params,'in_scale'));
                    texture.in_scale = params.in_scale;
                    texture.x0 = imresize(texture.x0,texture.in_scale,'bicubic');
                    texture.y0 = rand(size(texture.x0));
                end;
                if (isfield(params,'y'));texture.y0 = params.y0;end;
                if (isfield(params,'iter'));texture.iter = params.iter;end;
                if (isfield(params,'scales'));texture.scales = params.scales;end;
                if (isfield(params,'out_scale'));texture.out_scale = params.out_scale;end;
                
            end
        end
        
        function add_constraint(texture,constraint)
            n = numel(texture.constraints);
            texture.constraints{n+1} = constraint;
        end
        
        
        function run_variational_synthesis(texture)
            x_sz = [size(texture.x0,1),size(texture.x0,2)];
            texture.y = texture.y0;
            for s = numel(texture.scales):-1:1
                texture.x = imresize(texture.x0, texture.scales(s) *x_sz, 'bicubic');
                texture.y = imresize(texture.y, texture.out_scale *texture.scales(s) *x_sz,'bicubic');
                
%                 for k = 1:numel(texture.constraints)
%                     change_scale(texture.constraints{k});
%                 end
                
                for i = 1:texture.iter
                    
                    tic
                    enforce_all_constraints(texture);
                    disp(['Elapsed time for iteration ',num2str(i),' is ',num2str(toc)]);
                    
                    subplot(121);imshow(texture.x);title(['Exemplar at scale ',num2str(texture.scales(s)),', iteration ',num2str(i)]);
                    subplot(122);imshow(texture.y);title('Synthesis');drawnow;
                end
                
                
            end
            
        end
        
        function enforce_all_constraints(texture)
            for k = 1:numel(texture.constraints)
                constraint = texture.constraints{k};
                enforce_constraint(constraint);
            end
        end
        
        
        
    end
    
end
    
