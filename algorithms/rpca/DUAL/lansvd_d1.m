%this is a simplified version of lansvd that is specialized to compute
%the largest singular value and one principal singular vector of A
%
% Minming Chen, June 2009. Questions? v-minmch@microsoft.com ; 
% Arvind Ganesh (abalasu2@illinois.edu)
%
% Copyright: Perception and Decision Laboratory, University of Illinois, Urbana-Champaign
%            Microsoft Research Asia, Beijing


function [S1] = lansvd_d1(A)

[m, n] = size(A);
K = min([m n 10]);
K0 = K;

stop = 0;
nrestart = 0;
p_j1 = rand(m,1) - 0.5;

while ~stop
    beta_j = sqrt(p_j1'*p_j1);
    u_j = p_j1/beta_j;
    v_j1 = zeros(n,1);

    j = 1;
    U = u_j;
    V = [];
    alpha = zeros(K + 1, 1);
    beta = zeros(K + 1, 1);
    beta(1) = beta_j;
    while j <= K %Lanczos bi-diagonalizaion
        r_j = A'*u_j - beta_j*v_j1;
    
%         for i = 1:j - 1 %full-orthogonalization for r_j
%             v_i = V(:,i);
%             t = v_i'*r_j;
%             r_j = r_j - t*v_i;
%         end

%         %full-orthogonalization for r_j, in the matrix-vector product form
        if j > 1
            t = V'*r_j;
            r_j = r_j - V*t;
        end
        
        alpha_j = sqrt(r_j'*r_j);
        v_j = r_j / alpha_j;
        V = [V v_j];
        alpha(j) = alpha_j;

        p_j = A*v_j - alpha_j*u_j;
%         for i = 1:j %full-orhogonalization for p_j
%             u_i = U(:,i);
%             t = u_i'*p_j;
%             p_j = p_j - t*u_i;
%         end

%         %full-orthogonalization for p_j, in the matrix-vector product
%         form
        t = U'*p_j;
        p_j = p_j - U*t;
        
        beta_j = sqrt(p_j'*p_j);
        u_j = p_j / beta_j;
        U = [U u_j];
        beta(j + 1) = beta_j;
      
        v_j1 = v_j;
        j = j + 1;  
    end

    B = spdiags([alpha(1:K) beta(2:K+1)],[0 -1],K + 1,K);
    [U_B,S_B,V_B] = svd(full(B),0);%as K is usually small, this way of computing the SVD of B is harmless, although not optimal

    bot_U = U_B(end,1:K);%the bottom line of U_B reveals the error from the true left singular vector
    bot_V = V_B(end,1:K);%the bottom line of V_B reveals the error from the true right singular vector

    bot_U(1:3);
    bot_V(1:3);
    
    bnd_U = abs(bot_U);
    
    U1 = U*U_B(:,1);
    
    if bnd_U(1) <= 1E-8 %the principal singular value found
        stop = 1;
        
        S1 = S_B(1,1);
        V1 = V*V_B(:,1);
    else
%        p_j1 = p_j;%restart
%        K = min([m n 2*K])
        p_j1 = U1;%use the Ritz vector to restart; much faster than using p_j to restart
        K = min([m n K + K0]);
%        K = max(K + K0, ceil(1.5*K))
    
        nrestart = nrestart + 1;
    end
end 
