function [  Z, H, dnorm, AC, MIhat  ] = seminmf( X, k, varargin )
% Matrix sizes
% X: m x n
% Z: m x num_of_components
% H: num_of_components x num_of_components

% Process optional arguments
pnames = {'z0' 'h0' 'bUpdateH' 'maxiter' 'TolFun' 'bUpdateZ' 'verbose'};

% Do SVD initialisation of the init components

[z0, h0] = NNDSVD(abs(X), k, 0);
%rng(0);
%h0 = abs(rand(k , size(X, 2)));
%z0 = X * pinv(h0) + eps;


% z0 = rand(size(X, 1), k);


dflts  = {z0, h0, 1, 300,  1e-5, 1, 1};

[Z, H, bUpdateH, max_iter, tolfun, bUpdateZ, verbose] = ...
        internal.stats.parseArgs(pnames,dflts,varargin{:});

for i = 1:max_iter
    
    if bUpdateZ
        Z = X * pinv(H);
    end
    
    A = Z' * X;
    Ap = (abs(A)+A)./2;
    An = (abs(A)-A)./2;
    
    B = Z' * Z;
    Bp = (abs(B)+B)./2;
    Bn = (abs(B)-B)./2;
    
    if bUpdateH
        H = H .* sqrt((Ap + Bn * H) ./ (An + Bp * H + eps));
    end
      
    if mod(i, 10) == 0 || mod(i+1, 10) == 0 
        
        s = X - Z * H;
        dnorm = sqrt(sum(s(:).^2));
        % dnorm = norm(gX - Z * H, 'fro');
        
        if mod(i+1, 10) == 0
            dnorm0 = dnorm;
            continue
        end

        if mod(i, 100) == 0 && verbose
            display(sprintf('...Semi-NMF iteration #%d out of %d, error: %f\n', i, max_iter, dnorm));
        end

        if 0 && exist('dnorm0')
            assert(dnorm <= dnorm0, sprintf('Rec. error increasing! From %f to %f. (%d)', dnorm0, dnorm, k));
        end

        % Check for convergence
        if exist('dnorm0') && dnorm0-dnorm <= tolfun*max(1,dnorm0)
            if verbose
                display(sprintf('Stopped at %d: dnorm: %f, dnorm0: %f', i, dnorm, dnorm0));
            end
            break;
        end
     
    end
end

