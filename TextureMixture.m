classdef TextureMixture < handle
    
    properties (SetAccess = public)
        x1_in
        x2_in
        y0
        iter
        scales
        out_scale
        type
        constraints = {};
        band_weights; 
        x1
        x2
        y
    end
    
    
    
    methods
        function texture = TextureMixture(x1_in,x2_in,varargin)
            texture.x1_in = x1_in;
            texture.x2_in = x2_in;
            %default parameters
            texture.y0 =  rand(size(texture.x1_in));
            texture.iter = 8;
            texture.scales = [1,.5,.25];
            texture.out_scale = 1;
            texture.type = 'alternating';
            texture.band_weights = [.7,.3,.7,.6]; %hl,lh,hh,ll
            
            if numel(varargin)
                params = varargin{1};
                read_params(texture,params);
            end
        end
        
        function add_constraint(texture,constraint)
            n = numel(texture.constraints);
            texture.constraints{n+1} = constraint;
        end
        
        
        function run_variational_synthesis(texture)
            x_sz = [size(texture.x1_in,1),size(texture.x1_in,2)];
            texture.y = texture.y0;
            for s = numel(texture.scales):-1:1
                texture.x1 = imresize(texture.x1_in, texture.scales(s) *x_sz, 'bicubic');
                texture.x2 = imresize(texture.x2_in, texture.scales(s) *x_sz, 'bicubic');
                
                texture.y = imresize(texture.y, texture.out_scale *texture.scales(s) *x_sz,'bicubic');

                
                for i = 1:texture.iter
                    
%                     tic
                    enforce_all_constraints(texture);
%                     disp(['Elapsed time for iteration ',num2str(i),' is ',num2str(toc)]);
                    
                    subplot(131);imshow(texture.x1);title(['Exemplar at scale ',num2str(texture.scales(s)),', iteration ',num2str(i)]);
                    subplot(132);imshow(texture.x2);title(['Exemplar at scale ',num2str(texture.scales(s)),', iteration ',num2str(i)]);
                    subplot(133);imshow(texture.y);title('Synthesis');drawnow;
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
    
