% Serhat Selcuk Bucak, bucakser@msu.edu
function  [W_new, h, A, B] = inmf( V_new, W, h, A, B, rdim, beta, alpha, maxiter)
%%
%[W_new, h, A, B] = inmf( V_new, W, h, A, B, rdim, beta, alpha, maxiter)
%INPUTS:
%V_new: new data sample, a column vector (d x 1)
%W: Mixing Matrix -or matrix of basis vectors- (d x rdim)
%h: initialization for the new encoding vector, a column vector (rdim x 1),
%A: Matrix to store cummulative V*H' 
%B: Matrix to store cummulative H*H'
%rdim: factorization rank
%beta: weighting coefficient for contribution of initial samples (i.e. A=beta*A+alpha*V_new*h';B=beta*B+alpha*h*h')
%alpha: weighting coefficient for contribution of the new sample
%maxiter: maximum number of iterations
%OUTPUTS
%W_new: Updated mixing matrix 
%h: Updated encoding vector
%A, B: Updated A,B matrices

% Dimensions
vdim = size(W,1);
%Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W_new = W;   % Warm start: Use the W obtained from the previous step as initial matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Start iteration
iter = 0;
while 1
  
    if (iter==maxiter)
        break
    end
    
    iter = iter+1;    
    % Save old values
    %Wold = W_new;
    %Hold = H;
    
    % Compute new W and H 
    h = h.*((W_new)'*V_new)./((W_new)'*(W_new)*h + 1e-9); 
    W_new = W_new.*((beta*A+alpha*V_new*(h)')./(W_new*(beta*B+alpha*(h)*(h)') + 1e-9)); % Add 1e-9 to avoid 0 in the denom.

    if iter>1
        oldobj = newobj;
    end
    newobj = ((sum(sum((V_new-W_new*h).^2)))/(vdim));
    
    if (iter>1 && (oldobj-newobj)/newobj<1e-4)
        break;
    end
   
end

A=beta*A+alpha*V_new*h';
B=beta*B+alpha*h*h';
