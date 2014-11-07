function op = prox_nuclear( q, LARGESCALE )

%PROX_NUCLEAR    Nuclear norm.
%    OP = PROX_NUCLEAR( q ) implements the nonsmooth function
%        OP(X) = q * sum(svd(X)).
%    Q is optional; if omitted, Q=1 is assumed. But if Q is supplied, 
%    it must be a positive real scalar.
%
%    OP = PROX_NUCLEAR( q, LARGESCALE )
%       uses a Lanczos-based SVD if LARGESCALE == true,
%       otherwise it uses a dense matrix SVD
%
%    CALLS = PROX_NUCLEAR( 'reset' )
%       resets the internal counter and returns the number of function
%       calls
%
% This implementation uses a naive approach that does not exploit any
% a priori knowledge that X and G are low rank or sparse. Future
% implementations of TFOCS will be able to handle low-rank matrices 
% more effectively.
% Dual: proj_spectral.m
% See also prox_trace.m  and proj_spectral.m

if nargin == 1 && strcmpi(q,'reset')
    op = prox_nuclear_impl;
    return;
end

if nargin == 0,
	q = 1;
elseif ~isnumeric( q ) || ~isreal( q ) || numel( q ) ~= 1 || q <= 0,
	error( 'Argument must be positive.' );
end
if nargin < 2, LARGESCALE = []; end

% clear the persistent values:
prox_nuclear_impl();

op = @(varargin)prox_nuclear_impl( q, LARGESCALE, varargin{:} );

end % end of main function

function [ v, X ] = prox_nuclear_impl( q, LARGESCALE, X, t )
persistent oldRank
persistent nCalls
if nargin == 0, oldRank = []; v = nCalls; nCalls = []; return; end
if isempty(nCalls), nCalls = 0; end

ND = (size(X,2) == 1);
% ND = ~ismatrix(X);
if ND, % X is a vector, not a matrix, so reshape it 
    sx = size(X);
    X = reshape( X, prod(sx(1:end-1)), sx(end) );
end

if nargin > 3 && t > 0,
    
    if ~isempty(LARGESCALE) % inherited from parent
        largescale = LARGESCALE;
    else
        largescale = ( numel(X) > 100^2 ) && issparse(X);
    end
    tau = q*t;
    nCalls = nCalls + 1;
    
%     fprintf('ranks: ');
    if ~largescale
        [U,S,V] = svd( full(X), 'econ' );
    else
        % Guess which singular value will have value near tau:
        [M,N] = size(X);
        if isempty(oldRank), K = 10;
        else, K = oldRank + 2;
        end
        
        ok = false;
        opts = [];
        opts.tol = 1e-10; % the default in svds
        opt  = [];
        opt.eta = eps; % makes compute_int slow
%         opt.eta = 0;  % makes reorth slow
        opt.delta = 10*opt.eta;
        while ~ok
            K = min( [K,M,N] );
            if exist('lansvd','file')
                [U,S,V] = lansvd(X,K,'L',opt );
            else
                [U,S,V] = svds(X,K,'L',opts);
            end
            ok = (min(diag(S)) < tau) || ( K == min(M,N) );
%             fprintf('%d ',K );
            if ok, break; end
%             K = K + 5;
            K = 2*K;
            if K > 10
                opts.tol = 1e-6;
            end
            if K > 40
                opts.tol = 1e-4;
            end
            if K > 100
                opts.tol = 1e-1;
            end
            if K > min(M,N)/2
%                 disp('Computing explicit SVD');
                [U,S,V] = svd( full(X), 'econ' );
                ok = true;
            end
        end
        oldRank = length(find(diag(S) > tau));
    end
    s  = diag(S) - tau;
    tt = s > 0;
    s  = s(tt,:);
    
%     fprintf('\n')';
%     fprintf('rank is %d\n', length(tt) );
    
    % Check to make sure this doesn't break existing code...
    if isempty(s),
%         X(:) = 0;  % this line breaks the packSVD version
        X = tfocs_zeros(X);
    else
        X = U(:,tt) * bsxfun( @times, s, V(:,tt)' );
    end
else
    s = svd(full(X)); % could be expensive!
end

v = q * sum(s);
if ND, 
    X = reshape( X, sx ); 
end

end


% TFOCS v1.3 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2013 California Institute of Technology and CVX Research.
% See the file LICENSE for full license information.
