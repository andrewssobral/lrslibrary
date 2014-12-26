function   [A, D, alpha] = LowRankDictionaryUpdate_m(A_prev, D_prev, alpha_prev, lambda)
 
% min \|E\|_1+ \lambda*(\|D\|_f^2+\|\alpha\|_f^2) s.t. D*\lambda+E=X

% this function is to update the bases D and coefficients \alpha, over
% the constraint D*\alpha = A_prev, by one iteration of ADM

% Xianbiao Shu (xshu2@illinois.edu)
% Updated on Aug-10-2011 
% Copyright: Mitsubishi Electric Research Lab

    D = A_prev*alpha_prev'*inv(alpha_prev*alpha_prev'+lambda);
    alpha = inv(D'*D+ lambda)*D'*A_prev;
    
    
    A = D*alpha;
    
end
    