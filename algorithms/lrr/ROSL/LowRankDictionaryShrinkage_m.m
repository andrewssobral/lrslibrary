

function   [A, D, alpha] = LowRankDictionaryShrinkage_m(A_prev, D_prev, alpha_prev, lambda)
% min \|E\|_1+ \lambda\sum_{1<i<T}\|\alpha_i\|_2, s.t. \|D_i\|_2= 1 ; 

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
alpha_norm = zeros(K,1);

    for i =1:K
        
        % compute Error
        D(:,i) = 0;
        Error = A_prev-D*alpha;
        
        % compute D(:,i)
        D(:,i) = Error*alpha(i,:)';
        D_norm = norm(D(:,i));
        
        
        
        
        if  D_norm >0
            
            
            
            % Gram-Schmidt on D
            for j=1:i-1
                
                D(:,i) = D(:,i)-D(:,j)*(D(:,j)'*D(:,i));  
            end
            
            % normalization
%             D(:,i)= D(:,i)/D_norm; 
            D(:,i)= D(:,i)/norm(D(:,i));
            

            % compute alpha(i,:)
            alpha(i,:) = D(:,i)'*Error;
            
            
            
            % Gram-Schmidt on alpha
%             for j=1:i-1
%                 alpha_norm = norm(alpha(j,:));
%                 if(alpha_norm>0)
%                     alpha(i, :) = alpha(i, :)- alpha(i, :)*alpha(j, :)'*alpha(j, :)/alpha_norm^2;
%                 end
%             end

             
            
            
            % Magnitude thresholding
            alpha_norm(i) = norm(alpha(i,:)); 
            
            alpha(i,:) = alpha(i,:)*max(alpha_norm(i)-lambda, 0)/alpha_norm(i);
            alpha_norm(i) = max(alpha_norm(i)-lambda, 0);

             
              
        else
            alpha(i,:)=0;
            alpha_norm(i) = 0;
        end
        
        
    end
        
    % delete the zero bases;
        i = 1; 
        while(i <= K)
            if (alpha_norm(i) == 0)
                D(:, i) =[];
                alpha(i, :)= [];
                alpha_norm(i)= [];
                K = K-1;
            else
                i= i+1;
            end 
            
        end
        
%     [U, S, V]= svd(alpha, 'econ'); 
%     diagS = diag(S);
%     svp = length(find(diagS > lambda));
%     alpha =  diag(diagS(1:svp) - lambda) * V(:, 1:svp)'; 
%     D = D * U(:, 1:svp);
    
    
%     [U, S, V]= svd(alpha, 'econ'); 
%     diagS = diag(S);
%     svp = length(find(diagS > 0));
%     alpha =  diag(diagS(1:svp)-0) * V(:, 1:svp)'; 
%     D = D * U(:, 1:svp);


    
    A = D*alpha;
    
end
    