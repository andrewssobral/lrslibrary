

function   [A, D, alpha] = LowRankDictionarySparsify_m(A_prev, D_prev, alpha_prev, lambda)
% min \|E\|_1+ \|\alpha\|_1, s.t. \|D_i\|_2= 1 ; 

% this function is to update the bases D and coefficients \alpha, over
% the constraint D*\alpha = A_prev, by stepwise Gauss-Seidal

% Xianbiao Shu (xshu2@illinois.edu)
% Copyright: Mitsubishi Electric Research Lab
% Reference: X. Shu, F. Porikli, N. Ahujia "Robust Orthonormal Subspace Learning: Efficient Recovery of Corrupted Low-rank Matrices". CVPR 2014; 

D = D_prev; 
alpha = alpha_prev;



% D(:,1) = 0;
% Error = A_prev-D*alpha;

K = size(D_prev, 2);


    for i =1:K
        
        % compute Error
        D(:,i) = 0;
        Error = A_prev-D*alpha;
        
        % compute D(:,i)
        D(:,i) = Error*alpha(i,:)';
        D_norm = norm(D(:,i));
        
        
        
        
        if  D_norm >0
            
    
            % normalization
            D(:,i)= D(:,i)/D_norm; 
            
            % compute alpha(i,:)
            alpha(i,:) = D(:,i)'*Error;
            
            
            % Soft thresholding
            
            alpha(i,:) = sign(alpha(i,:)).*max(abs(alpha(i,:))-lambda, 0);
        else
            alpha(i,:)=0;
        end
        
        
    end
        

      
    
    


    
    A = D*alpha;
    
end
    