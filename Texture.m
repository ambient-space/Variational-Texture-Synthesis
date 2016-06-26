classdef Texture < handle
    
    properties (SetAccess = public)
        x0
        y0
        iter
        scales
        out_scale
        type
        constraints = {};
        x
        y
        
        %for graceful timeout
        start_time
        algorithm_time
        time_remaining
        current_time_remaining
        time_per_scale
    end
    
    
    
    methods
        function texture = Texture(x_in,varargin)
            texture.x0 = x_in;
            %default parameters
            texture.y0 =  rand(size(texture.x0));
            texture.iter = 8;
            texture.scales = [1,.5,.25];
            texture.out_scale = 1;
            texture.type = 'alternating';
            
            if numel(varargin)
                params = varargin{1};
                read_params(texture,params);
            end
            
            tmp = texture.scales .^2;
            texture.time_per_scale  = tmp/sum(tmp(:));
            
        end
        
        function add_constraint(texture,constraint)
            n = numel(texture.constraints);
            texture.constraints{n+1} = constraint;
        end
        
        
        function run_variational_synthesis(texture)
            x_sz = [size(texture.x0,1),size(texture.x0,2)];
            texture.y = texture.y0;
            texture.start_time = tic;
            
            
            for s = numel(texture.scales):-1:1
%                 texture.current_time_remaining = texture.time_per_scale * (texture.algorithm_time - toc(texture.start_time));
                
                texture.x = imresize(texture.x0, texture.scales(s) *x_sz, 'bicubic');
                texture.y = imresize(texture.y, texture.out_scale *texture.scales(s) *x_sz,'bicubic');

                
                for i = 1:texture.iter
                    
%                     tic
                    enforce_all_constraints(texture);
%                     disp(['Elapsed time for iteration ',num2str(i),' is ',num2str(toc)]);
                    
                    subplot(121);imshow(texture.x);title(['Exemplar at scale ',num2str(texture.scales(s)),', iteration ',num2str(i)]);
                    subplot(122);imshow(texture.y);title('Synthesis');drawnow;
                end
                
                
            end
            
            tmp_sz = size(texture.y0);
            texture.y = imresize(texture.y,tmp_sz(1:2));
            
        end
        
        function enforce_all_constraints(texture)
            for k = 1:numel(texture.constraints)
                constraint = texture.constraints{k};
                enforce_constraint(constraint);
            end
        end
        
        
    end
    
end
    
