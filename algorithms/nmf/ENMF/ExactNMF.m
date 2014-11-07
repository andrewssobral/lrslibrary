function [V,W] = ExactNMF(S, r, max_attempts)
% Tries to compute an exact rank-r nonnegative factorization of matrix S
% (each line is a different attempt, and displays the decreasing approximation 
% errors - first number in brackets is best error seen so far)
% If not found, returns the best solution found so far. 
% 
% This is the method used to compute the factorizations in the paper: 
% N. Gillis and F. Glineur, "On the Geometric Interpretation of the 
% Nonnegative Rank", Linear Algebra and its Applications 437 (11), 
% pp. 2685-2712, 2012.
%
% [V,W] = ExactNMF(S, r, max_attempts)
%
% Input.
%   S              : (m x n) nonnegative matrix to factorize
%   r              : factorization rank
%   max_attempts   : maximum numberof attemps to find an exact NMF
%
% Output.
%   (V,W)    : nonnegative matrices (m x r) and (r x n) ||VW-M|| is 
%              smallest among all max_attempts initializations (except 
%              if an exact factorization earlier), hopefully zero...

tic;
smallest=inf;
if nargin<3
   max_attempts = 100;
end
for i=1:max_attempts  
   V=rand(size(S,1),r);W=rand(r,size(S,2));
   [V,W]=HALS(S,V,W,0); % apply optimal scaling   
   olderr=norm(S-V*W,'fro');
   fprintf('\n.#%4d (%6.3f) - ',i, smallest);
   while 1
   fprintf('%6.3f ',olderr);
   [V,W]=HALS(S,V,W,100);
   err=norm(S-V*W,'fro');
   if err<1e-6
      [V,W]=HALS(S,V,W,1000); % refine solution
      fprintf('\nEXACT factorization found after %d attempts (error=%g)\n', i, err);
      return;
   elseif err>olderr/2
      if err<smallest
         smallest=err;
         smallestV=V; smallestW=W;
      end
      break;
   end
   olderr=err;
   end
end
[V,W]=HALS(S,smallestV,smallestW,1000);
fprintf('\nNo exact factorization found in %d attempts (%d secs) - best error is %g.\n', max_attempts, toc, norm(S-V*W,'fro'));
disp('Please, try again...'); 