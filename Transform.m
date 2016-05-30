classdef Transform
    properties (Constant)
        %If you can compile MEX, im2col,col2im perform
        %more than 50x faster for large arrays
        use_mex = 1;
    end
    
    methods (Static)
        function X = im2col_sample(x,blocksize,ids)
            if Transform.use_mex
                X = im2col_sample_mex(x,ids,blocksize,[size(x,1),size(x,2),size(x,3)]);
            else
                N=size(x,3);
                sz=[size(x,1),size(x,2)];
                X=zeros([N*blocksize^2,numel(ids)],'double');
                r=0:blocksize-1;
                for k=1:numel(ids)
                    [i,j]=ind2sub(sz,ids(k));
                    X(:,k)=reshape(x(r+i,r+j,:),[size(X,1),1]);
                end
            end
        end
        
        function x = col2im_sample(X,x,blocksize,ids_x)
            if Transform.use_mex
                sz = [size(x,1),size(x,2),size(x,3)];
                tmp_x = col2im_sample_mex(X,ids_x,blocksize,sz);
            else
                N=size(X,1)/blocksize^2;
                if N~=fix(N)
                    error('incorrect dims');
                end
                x_sz = size(x);
                tmp_x = zeros(x_sz);
                x_sz = x_sz(1:2);
                weightmatrix = zeros(x_sz(1:2));
                for k=1:size(X,2)
                    [i,j] = ind2sub(x_sz,ids_x(k));
                    temp = X(:,k);
                    weightmatrix(i:i+blocksize-1,j:j+blocksize-1)=weightmatrix(i:i+blocksize-1,j:j+blocksize-1)+1;
                    tmp_x(i:i+blocksize-1,j:j+blocksize-1,:) = tmp_x(i:i+blocksize-1,j:j+blocksize-1,:)+reshape(temp,blocksize,blocksize,N);
                end
                
                weightmatrix(weightmatrix~=0)=1./weightmatrix(weightmatrix~=0);
                
                for i=1:size(tmp_x,3)
                    tmp_x(:,:,i) = tmp_x(:,:,i).*weightmatrix;
                end
            end
            x(tmp_x~=0) = tmp_x(tmp_x~=0); %image may not be covered in sample
        end
        
        
        
        function [x_pca,base] = pca(x,varargin)
            x_vec = reshape(x,[],size(x,3))';
            %varargin optional param specifying basis
            if numel(varargin)
                base = varargin{1};
            else
                [base, ~] = eig(x_vec * x_vec');
            end
            
            % change to pca basis
            x_pca = reshape((base' * x_vec)',size(x));
        end
        
        function x = ipca(x_pca,base)
            x_vec = reshape(x_pca,[],size(x_pca,3))';
            x = reshape((base*x_vec)',size(x_pca));
        end
        
        function x_udwt = haar_udwt_2D(x,num_levels)
            %filter RGB seperately
            if size(x,3) > 1
                x_udwt = [];
                for c = 1:size(x,3)
                    x_udwt = cat(3,x_udwt,Transform.haar_udwt_2D(x(:,:,c),num_levels));
                end
            else
                x_sz = size(x);
                x_udwt = zeros([x_sz,num_levels*3+1]);
                curr_dec = 1;
                for l = 1:num_levels
                    band_num = (l-1)*3 +1;
                    if l == 1
                        x_ll = x;
                    end
                    
                    x1 =  circshift(x_ll,[0,curr_dec]);
                    x2 =  circshift(x_ll,[curr_dec,0]);
                    x3 =  circshift(x_ll,[curr_dec,curr_dec]);
                    
                    x_hl = x_ll + x1 - x2 -x3;
                    x_lh = x_ll - x1 + x2 -x3;
                    x_hh = x_ll - x1 - x2 + x3;
                    
                    x_udwt(:,:,band_num) = 1/4*x_hl;
                    x_udwt(:,:,band_num+1) = 1/4*x_lh;
                    x_udwt(:,:,band_num+2) = 1/4*x_hh;
                    
                    x_ll = 1/4*(x_ll + x1 + x2 + x3);
                    
                    
                    if l == num_levels
                        x_udwt(:,:,end) = x_ll;
                    end
                    curr_dec = curr_dec*2;
                end
                
                
            end
        end
        
        
        function x = haar_iudwt_2D(x_sz,x_udwt)
            
            if numel(x_sz) == 3
                C = x_sz(3);
                num_levels = (size(x_udwt,3)-C)/(3*C);
                x = [];
                for c = 1:C
                    band_ids = (c-1)*(3*num_levels+1)+1 : c*(3*num_levels+1);
                    x = cat(3,x,Transform.haar_iudwt_2D(x_sz(1:2),x_udwt(:,:,band_ids)));
                end
            else
                n = (size(x_udwt,3)-1)/3;
                x_ll = x_udwt(:,:,end);
                for l = n:-1:1
                    band_num = (l-1)*3+1;
                    curr_dec = -2^(l-1);
                    
                    x_band = x_udwt(:,:,band_num);
                    x1 =  circshift(x_band,[0,curr_dec]);
                    x2 =  circshift(x_band,[curr_dec,0]);
                    x3 =  circshift(x_band,[curr_dec,curr_dec]);
                    x_hl = x_band +x1 - x2 -x3;
                    
                    x_band = x_udwt(:,:,band_num+1);
                    x1 =  circshift(x_band,[0,curr_dec]);
                    x2 =  circshift(x_band,[curr_dec,0]);
                    x3 =  circshift(x_band,[curr_dec,curr_dec]);
                    x_lh = x_band -x1 +x2 -x3;
                    
                    x_band = x_udwt(:,:,band_num+2);
                    x1 =  circshift(x_band,[0,curr_dec]);
                    x2 =  circshift(x_band,[curr_dec,0]);
                    x3 =  circshift(x_band,[curr_dec,curr_dec]);
                    x_hh = x_band - x1 - x2 + x3;
                    
                    x_band = x_ll;
                    x1 =  circshift(x_band,[0,curr_dec]);
                    x2 =  circshift(x_band,[curr_dec,0]);
                    x3 =  circshift(x_band,[curr_dec,curr_dec]);
                    x_ll = x_ll + x1 + x2 + x3;
                    
                    x_ll = 1/4*(x_ll + x_hl + x_lh + x_hh);
                end
                x = x_ll;
            end
        end
        
        
        
        
    end
    
    
end

