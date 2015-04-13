function [A,E,Y] = singular_value_rpca( D, lambda, tau, delta, svdMethod, A0)

%   Solves the Robust PCA relaxation
%
%      min  \tau ( |A|_* + \lambda |E|_1 ) + 1/2 |(A,E)|_F^2
%      subj  A+E = D
%
%   by iterative thresholding.
%
%   Inputs:
%      D      -- the data matrix, m x n.
%      lambda -- relative weight of sparsity of error
%
%   Optional:
%      tau    -- magnitude of L2 relaxation of the pure robust PCA SDP,
%                 higher value is desirable.
%      delta  -- stepsize, should be in (0,1).
%      svdMethod -- SVD routine to be used in each iteration, must be one
%                    of 'svdlibc', 'propack', or 'svds'. If not any of
%                    these, default option is MATLAB's svd command. (May
%                    require additional library dependencies if custom routine is used.)
%      A0 -- true low-rank solution, if known, to enable better display of
%            progress in each iteration.
%
%   Outputs:
%      A      -- estimate of the low-rank generating matrix
%      E      -- estimate of the error or corruption
%
%   Winter '08, John Wright, Shankar Rao. Questions? jnwright@uiuc.edu
%
%   Copyright: Perception and Decision Laboratory
%			University of Illinois, Urbana-Champaign 

VERBOSE                         = 2;
EPSILON_PRIMAL                  = 5e-4;

if nargin < 5,  svdMethod = 'svd'; end

if nargin < 4,  delta = 0.9; end;

if nargin < 3, tau = 1e4; end;


MAX_ITER                        = 25000;
DISPLAY_EVERY                   = 100;

[m,n] = size(D);

Y = zeros(m,n);  % Lagrange multiplier
A = zeros(m,n);  % Structure
E = zeros(m,n);  % Error

rankA = 0;

iter      = 0;
converged = false;

while ~converged
    iter = iter + 1;
    
    switch lower(svdMethod)
        case 'svdlibc'
            [U,diagS,V] = svdlibc(Y, rankA+1);
        case 'propack'
            [U,S,V] = lansvd(Y,rankA+1,'L');
            diagS = diag(S);
        case 'svds'
            [U,S,V] = svds(Y, rankA+1, 'L');
            diagS = diag(S);
        otherwise            
            [U,S,V] = svd(Y,0);            
            diagS = diag(S);            
    end
    
    
    A = U * diag(pos(diagS-tau)) * V';
    E = sign(Y) .* pos( abs(Y) - lambda*tau );
    M = D - A - E;
    
    rankA  = sum(diagS>tau);
    cardE = sum(sum(double(abs(E)>0)));
    
    Y = Y + delta * M;
    
%     if VERBOSE > 1 && mod(iter, DISPLAY_EVERY)==0 && nargin>=6,
%         disp(['    Iteration '    num2str(iter)               ...
%             ' |A|_F '         num2str(norm(A,'fro'))      ...
%             ' rank(A) '         num2str(rankA)              ...
%             ' |E|_F '         num2str(norm(E,'fro'))      ...
%             ' |E|_0 '         num2str(cardE) ...
%             ' |D-A-E|_F ' num2str(norm(M,'fro')) ...
%             ' |A-A0|_F / |A0|_F ' num2str(norm(A-A0,'fro')/norm(A0,'fro')) ...
%             ' |D-A-E|_1,inf ' num2str(max(max(abs(M)))) ]);
%     elseif VERBOSE > 0 && mod(iter, DISPLAY_EVERY)==0,
%         disp(['    Iteration '    num2str(iter)               ...
%             ' rank(A) '         num2str(rankA) ...
%             ' ||E||_0 ' num2str(cardE) ]);
%     end
    
    if ( norm(D-A-E,'fro')/norm(D,'fro') < EPSILON_PRIMAL || iter >= MAX_ITER )
        converged = true;
    end
    
%     if ( iter >= MAX_ITER )
%         disp('Maximum number of iterations reached.') ;
%     end
      
end