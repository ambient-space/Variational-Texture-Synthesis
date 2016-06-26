classdef Transport
    properties (Constant)
        use_gpu = 0;
    end
    
    methods (Static)
        
        function [Y,NN_ids] = nn_search(Y,X)
%             fprintf('Computing nearest neighbors... ');
%             t0 = tic;
            %lower bound on OT problem
            if Transport.use_gpu
            	[Y,NN_ids] = Transport.nn_search_gpu(Y,X);
            else
            	[Y,NN_ids] = Transport.nn_search_cpu(Y,X);
            end
%             t = toc(t0);
%             fprintf(['Elapsed time is: ',num2str(t),'\n']);
        end
        
        function [Y,NN_ids] = nn_search_cpu(Y,X)
            X_l2 = (-2)*[ X ;ones( 1 , size(X,2) ); -1/2 * sum( X.^2 )  ];
            
            Y_l2 = [ Y ; -1/2 * sum( Y.^2 ) ; ones( 1 , size(Y,2) ) ];
            
            X_l2 = single(X_l2);
            Y_l2 = single(Y_l2);
            
            B = 100;
            
            N = ceil(size(Y_l2,2)/B);
            Y_l2 = [Y_l2,zeros(size(Y_l2,1),N*B-size(Y_l2,2))];%pad Y_l2
            NN_ids = zeros(B,N);
            for k =1:N
                r = (k-1)*B+1:k*B;
                D = X_l2'*Y_l2(:,r);
                [~,tmp] = min(D,[],1);
                NN_ids(:,k) = tmp(:);
            end
            NN_ids = NN_ids(1:size(Y,2));
            
            Y = X(:,NN_ids);
        end
        
        
        function [Y,NN_ids] = nn_search_gpu(Y,X)
            X_l2 = (-2)*[ X ;ones( 1 , size(X,2) ); -1/2 * sum( X.^2 )  ];
            
            Y_l2 = [ Y ; -1/2 * sum( Y.^2 ) ; ones( 1 , size(Y,2) ) ];
            
            X_l2 = gpuArray(single(X_l2));
            Y_l2 = gpuArray(single(Y_l2));
            
            B = 1500;
            N = ceil(size(Y_l2,2)/B);
            Y_l2 = [Y_l2,zeros(size(Y_l2,1),N*B-size(Y_l2,2))];%pad Y_l2
            NN_ids = gpuArray(single(zeros(B,N)));
            for k =1:N
                r = (k-1)*B+1:k*B;
                D = X_l2'*Y_l2(:,r);
                
                [~,tmp] = min(D,[],1);
                NN_ids(:,k) = tmp(:);
            end
            NN_ids = NN_ids(1:size(Y,2));
            
            NN_ids = gather(NN_ids);
            
            Y = X(:,NN_ids);

        end
        
        function y = optimal_transport_1D(y,x)
            [~,y_ids] = sort(y(:));
            [x_sort,~] = sort(x(:));
            x_sort = imresize(x_sort,[numel(y),1],'bilinear'); %bilinear interp
            y(y_ids) = x_sort;
        end
        
        
        function Y = sparsecode(X,D,sparsity)
            %This mex function needs to be compiled to work
            Y = omp_chol(D'*X,D'*D,sparsity);
        end
        
        
        
        
%         function y = transfer_autocorr(y,x)
%             Na = 5;
%             la = floor((Na-1)/2);
%             ace = NaN * ones(Na,Na);
%             ch = x;
%             [Nly, Nlx] = size(ch);
%             Sch = min(Nlx, Nly);
%             le = min(Sch/2-1,la);
%             cx = Nlx/2+1;  %Assumes Nlx even
%             cy = Nly/2+1;
%             ac = fftshift(real(ifft2(abs(fft2(ch)).^2)))/prod(size(ch));
%             ac = ac(cy-le:cy+le,cx-le:cx+le);
%             ace(la-le+1:la+le+1,la-le+1:la+le+1) = ac;
%             
%             [y, ~] = modacor22(y, ace,1);
%         end
        
        
        
    end
    
    
end

