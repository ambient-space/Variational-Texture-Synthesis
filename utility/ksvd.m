function [D,Y] = ksvd(X,dictsize,sparsity)
%initialize dictionary as random signals

p = randperm(size(X,2));
D = X(:,p(1:dictsize));

iter = 15; %fixed No. of iterations
for atom = 1 : dictsize
    D(:,atom) = D(:,atom)/(sqrt(sum(D(:,atom).^2)));
end

Y = ompmex_qt(D'*X,D'*D,sparsity);

for i = 1:iter
    
%     disp(['KSVD iteration: ',num2str(i)]);
    
    Y = ompmex_qt(D'*X,D'*D,sparsity);
    
    for atom = 1 : dictsize
        
        Y_atom = Y(atom,:)';
        
        %extract indices of signals that use this atom
        signal_ids = find(Y_atom);
        
        if signal_ids > 1
            Y_atom = Y_atom(signal_ids);
            
            % Approx. SVD of patches using this atom. The current atom + a gradient
            D(:,atom) = D(:,atom)*(Y_atom' * Y_atom) + ( X(:,signal_ids)*Y_atom - ( D*Y(:,signal_ids) )*Y_atom );
            
            %normalize atom
            D(:,atom) = D(:,atom)/(sqrt(sum(D(:,atom).^2)));
        end
        
          
    end
    
    D = cleardict(D,Y,X);
    
end




end



function D = cleardict(D,Gamma,X)

use_thresh = 4;  % at least this number of samples must use the atom to be kept
muthresh = .98; % promotes diverse atoms
dictsize = size(D,2);
unused_sigs = 1:size(X,2);
% compute error in blocks to conserve memory
err = zeros(1,size(X,2));
blocks = [1:3000:size(X,2) size(X,2)+1];
for i = 1:length(blocks)-1
  err(blocks(i):blocks(i+1)-1) = sum((X(:,blocks(i):blocks(i+1)-1)-D*Gamma(:,blocks(i):blocks(i+1)-1)).^2);
end

usecount = sum(abs(Gamma)>1e-7, 2);

for j = 1:dictsize
  
  % compute G(:,j)
  Gj = D'*D(:,j);
  Gj(j) = 0;
  
  % replace atom
  if ( (max(Gj.^2)>muthresh^2 || usecount(j)<use_thresh)  )
    [y,i] = max(err(unused_sigs));
    D(:,j) = X(:,unused_sigs(i)) / norm(X(:,unused_sigs(i)));
    unused_sigs = unused_sigs([1:i-1,i+1:end]);
  end
end

end
