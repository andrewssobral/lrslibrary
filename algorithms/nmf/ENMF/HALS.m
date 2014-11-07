% HALS algorithm for NMF
%
% See, e.g., N. Gillis and F. Glineur, "Accelerated Multiplicative Updates 
% and Hierarchical ALS Algorithms for Nonnegative Matrix Factorization", 
% Neural Computation 24 (4), pp. 1085-1105, 2012.
%
% [V,W] = HALS(M,Vinit,Winit,maxiter)
%
% Input.
%   M              : (m x n) matrix to factorize
%   (Vinit, Winit) : initial matrices of dimensions (m x r) and (r x n)
%   maxiter        : number of iterations
%
% Output.
%   (V,W)    : nonnegative matrices (m x r) and (r x n) s.t. VW approximates M


function [V,W] = HALS(M,Vinit,Winit,maxiter)

[m,n] = size(M);
[m,r] = size(Vinit);

% Initialization
V=Vinit; W=Winit;

% Scaling
A = V*W; % 2mnr operations
alpha = sum(sum(A.*M))/sum(sum( (V'*V).*(W*W')  )); % O(mn) operations
V = V*alpha; % mr operations

epsilon = 1e-16;

for j = 1 : maxiter
    A = (M*W'); % 2mnr operations
    B = (W*W'); % 2mr^2 operations
    for k = 1 : r
            % r times 2mr operations + O(m)
            V(:,k) = max(A(:,k)-V(:,[1:k-1 k+1:end])*B([1:k-1 k+1:end],k),epsilon) / (B(k,k)+1e-16);
    end

    A = (V'*M); % 2mnr operations
    B = (V'*V); % 2nr^2 operations
    for k = 1 : r
        % r times 2nr operations + O(n)
        W(k,:) = max(A(k,:)- B(k,[1:k-1 k+1:end]) * W([1:k-1 k+1:end],:),epsilon) /(B(k,k)+1e-16);
    end
end
% Set small values to zero
V(V==1e-16) = 0;  W(W==1e-16) = 0; 

% TOTAL : 4mnr + 4(m+n)r^2 + O(r(m+n)) operations