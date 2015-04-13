function M = fixedrankfactory_3factors_preconditioned(m, n, k)
% Manifold of m-by-n matrices of rank k with polar quotient geometry.
%
% function M = fixedrankLSRquotientfactory(m, n, k)
%
% A point X on the manifold is represented as a structure with three
% fields: L, S and R. The matrices L (mxk) and R (nxk) are orthonormal,
% while the matrix S (kxk) is a full rank matrix
% matrix.
%
% Tangent vectors are represented as a structure with three fields: L, S
% and R.

% This file is part of Manopt: www.manopt.org.
% Original author: Bamdev Mishra, Dec. 30, 2012.
% Contributors:
    
    % Change log:
    
    M.name = @() sprintf('LSR'' quotient manifold of %dx%d matrices of rank %d', m, n, k);
    
    M.dim = @() (m+n-k)*k;
    
    % Some precomputations at the point X to be used in the inner product (and
    % pretty much everywhere else).
    function X = prepare(X)
        if ~all(isfield(X,{'StS','SSt','invStS','invSSt'}) == 1)
            X.SSt = X.S*X.S';
            X.StS = X.S'*X.S;
            X.invSSt = eye(size(X.S, 2))/X.SSt;
            X.invStS = eye(size(X.S, 2))/X.StS;
        end
    end
    
    % Choice of the metric on the orthnormal space is the low-rank matrix completio cost function.
    M.inner = @iproduct;
    function ip = iproduct(X, eta, zeta)
        X = prepare(X);
        
        ip = trace(X.SSt*(eta.L'*zeta.L)) + trace(X.StS*(eta.R'*zeta.R)) ...
            + trace(eta.S'*zeta.S);
        
    end
    M.norm = @(X, eta) sqrt(M.inner(X, eta, eta));
    
    M.dist = @(x, y) error('fixedrankLSRquotientfactory.dist not implemented yet.');
    
    M.typicaldist = @() 10*k;
    
    skew = @(X) .5*(X-X');
    symm = @(X) .5*(X+X');
    
    M.egrad2rgrad = @egrad2rgrad;
    function rgrad = egrad2rgrad(X, egrad)
        X = prepare(X);
        
        SSL = X.SSt;
        ASL = 2*symm(SSL*(egrad.S*X.S'));
        
        SSR = X.StS;
        ASR = 2*symm(SSR*(egrad.S'*X.S));
        
        %         BL1 = lyap(SSL, -ASL);
        %         BR1 = lyap(SSR, -ASR);
        [BL, BR] = tangent_space_lyap(X.S, ASL, ASR);
        
        rgrad.L = (egrad.L - X.L*BL)*X.invSSt;
        rgrad.R = (egrad.R - X.R*BR)*X.invStS;
        
        
        rgrad.S = egrad.S;
        
        %                         norm(skew(X.SSt*(rgrad.L'*X.L) + rgrad.S*X.S'), 'fro')
        %                         norm(skew(X.StS*(rgrad.R'*X.R) - X.S'*rgrad.S), 'fro')
        
    end
    
    
    
    M.ehess2rhess = @ehess2rhess;
    function Hess = ehess2rhess(X, egrad, ehess, eta)
        X = prepare(X);
        
        % Riemannian gradient
        SSL = X.SSt;
        ASL = 2*symm(SSL*(egrad.S*X.S'));
        SSR = X.StS;
        ASR = 2*symm(SSR*(egrad.S'*X.S));
        [BL, BR] = tangent_space_lyap(X.S, ASL, ASR);
        
        rgrad.L = (egrad.L - X.L*BL)*X.invSSt;
        rgrad.R = (egrad.R - X.R*BR)*X.invStS;
        rgrad.S = egrad.S;
        
        % Directional derivative of the Riemannian gradient
        ASLdot = 2*symm((2*symm(X.S*eta.S')*(egrad.S*X.S')) + X.SSt*(ehess.S*X.S' + egrad.S*eta.S')) - 4*symm(symm(eta.S*X.S')*BL);
        ASRdot = 2*symm((2*symm(X.S'*eta.S)*(egrad.S'*X.S)) + X.StS*(ehess.S'*X.S + egrad.S'*eta.S)) - 4*symm(symm(eta.S'*X.S)*BR);
        
        %         SSLdot = X.SSt;
        %         SSRdot = X.StS;
        %         BLdot = lyap(SSLdot, -ASLdot);
        %         BRdot = lyap(SSRdot, -ASRdot);
        
        [BLdot, BRdot] = tangent_space_lyap(X.S, ASLdot, ASRdot);
        
        Hess.L = (ehess.L - eta.L*BL - X.L*BLdot - 2*rgrad.L*symm(eta.S*X.S'))*X.invSSt;
        Hess.R = (ehess.R - eta.R*BR - X.R*BRdot - 2*rgrad.R*symm(eta.S'*X.S))*X.invStS;
        Hess.S = ehess.S;
        
        
        
        % BM comments: Till this, everything seems correct.
        % We still need a correction factor for the non-constant metric
        % The correction factor owes itself to the Koszul formula...
        % This is the Riemannian connection in the Euclidean space with the
        % scaled metric.
        Hess.L = Hess.L + (eta.L*symm(rgrad.S*X.S') + rgrad.L*symm(eta.S*X.S'))*X.invSSt;
        Hess.R = Hess.R + (eta.R*symm(rgrad.S'*X.S) + rgrad.R*symm(eta.S'*X.S))*X.invStS;
        Hess.S = Hess.S - symm(rgrad.L'*eta.L)*X.S - X.S*symm(rgrad.R'*eta.R);
        
        % The Riemannian connection on the quotient space is the
        % projection on the tangent space of the total space and then onto the horizontal
        % space. This is accomplished by the following operation.
        Hess = M.proj(X, Hess);
        
        %         norm(skew(X.SSt*(Hess.L'*X.L) + Hess.S*X.S'))
        %         norm(skew(X.StS*(Hess.R'*X.R) - X.S'*Hess.S))
        
    end
    
    
    
    
    M.proj = @projection;
    function etaproj = projection(X, eta)
        X = prepare(X);
        
        
        % First, projection onto the tangent space of the total sapce
        SSL = X.SSt;
        ASL = 2*symm(X.SSt*(X.L'*eta.L)*X.SSt);
        BL = lyap(SSL, -ASL);
        eta.L = eta.L - X.L*BL*X.invSSt;
        
        SSR = X.StS;
        ASR = 2*symm(X.StS*(X.R'*eta.R)*X.StS);
        BR = lyap(SSR, -ASR);
        eta.R = eta.R - X.R*BR*X.invStS;
        
        % Project onto the horizontal space
        PU = skew((X.L'*eta.L)*X.SSt) + skew(X.S*eta.S');
        PV = skew((X.R'*eta.R)*X.StS)  + skew(X.S'*eta.S);
        [Omega1, Omega2] = coupled_lyap(X.S, PU, PV);
        %         norm(2*skew(Omega1*X.SSt) - PU -(X.S*Omega2*X.S'),'fro' )
        %         norm(2*skew(Omega2*X.StS) - PV -(X.S'*Omega1*X.S),'fro' )
        %
        
        etaproj.L = eta.L - (X.L*Omega1);
        etaproj.S = eta.S - (X.S*Omega2 - Omega1*X.S) ;
        etaproj.R = eta.R - (X.R*Omega2);
        
        %                                 norm(skew(X.SSt*(etaproj.L'*X.L) + etaproj.S*X.S'))
        %                                 norm(skew(X.StS*(etaproj.R'*X.R) - X.S'*etaproj.S))
        %
        %                                  norm(skew(X.SSt*(etaproj.L'*X.L) - X.S*etaproj.S'))
        %                                 norm(skew(X.StS*(etaproj.R'*X.R) + etaproj.S'*X.S))
        
    end
    
    
    M.tangent = M.proj;
    M.tangent2ambient = @(X, eta) eta;
    
    M.retr = @retraction;
    function Y = retraction(X, eta, t)
        if nargin < 3
            t = 1.0;
        end
        
        Y.S = (X.S + t*eta.S);
        Y.L = uf((X.L + t*eta.L));
        Y.R = uf((X.R + t*eta.R));
        
        Y = prepare(Y);
    end
    
    M.exp = @exponential;
    function Y = exponential(X, eta, t)
        if nargin < 3
            t = 1.0;
        end
        Y = retraction(X, eta, t);
        warning('manopt:fixedrankLSRquotientfactory:exp', ...
            ['Exponential for fixed rank ' ...
            'manifold not implemented yet. Used retraction instead.']);
    end
    
    M.hash = @(X) ['z' hashmd5([X.L(:) ; X.S(:) ; X.R(:)])];
    
    M.rand = @random;
    % Factors L and R live on Stiefel manifolds, hence we will reuse
    % their random generator.
    stiefelm = stiefelfactory(m, k);
    stiefeln = stiefelfactory(n, k);
    function X = random()
        X.L = stiefelm.rand();
        X.R = stiefeln.rand();
        X.S = diag(1+rand(k, 1));
        
        X = prepare(X);
    end
    
    M.randvec = @randomvec;
    function eta = randomvec(X)
        % A random vector on the horizontal space
        eta.L = randn(m, k);
        eta.R = randn(n, k);
        eta.S = randn(k, k);
        eta = projection(X, eta);
        nrm = M.norm(X, eta);
        eta.L = eta.L / nrm;
        eta.R = eta.R / nrm;
        eta.S = eta.S / nrm;
    end
    
    M.lincomb = @lincomb;
    
    M.zerovec = @(X) struct('L', zeros(m, k), 'S', zeros(k, k), ...
        'R', zeros(n, k));
    
    M.transp = @(x1, x2, d) projection(x2, d);
    
    % vec and mat are not isometries, because of the unusual inner metric.
    M.vec = @(X, U) [U.L(:) ; U.S(:); U.R(:)];
    M.mat = @(X, u) struct('L', reshape(u(1:(m*k)), m, k), ...
        'S', reshape(u((m*k+1): m*k + k*k), k, k), ...
        'R', reshape(u((m*k+ k*k + 1):end), n, k));
    M.vecmatareisometries = @() false;
    
end

% Linear combination of tangent vectors
function d = lincomb(x, a1, d1, a2, d2) %#ok<INLSL>
    
    if nargin == 3
        d.L = a1*d1.L;
        d.R = a1*d1.R;
        d.S = a1*d1.S;
    elseif nargin == 5
        d.L = a1*d1.L + a2*d2.L;
        d.R = a1*d1.R + a2*d2.R;
        d.S = a1*d1.S + a2*d2.S;
    else
        error('Bad use of fixedrankLSRquotientfactory.lincomb.');
    end
    
end

function A = uf(A)
    [L, unused, R] = svd(A, 0); %#ok
    A = L*R';
end

function[BU, BV] = tangent_space_lyap(R, E, F)
    % We intent to solve     RR^T  BU + BU RR^T  = E
    %                        R^T R BV + BV R^T R = F
    %
    % This can be solved using two calls to the Matlab lyap.
    % However, we can still have a more efficient implementations as shown
    % below...
    
    [U, Sigma, V] = svd(R);
    E_mod = U'*E*U;
    F_mod = V'*F*V;
    b1 = E_mod(:);
    b2 = F_mod(:);
    
    r = size(Sigma, 1);
    sig = diag(Sigma); % all the singular values in a vector
    sig1 = sig*ones(1, r); % columns repeat
    sig1t = sig1'; % rows repeat
    s1 = sig1(:);
    s2 = sig1t(:);
    
    % The block elements
    a =  s1.^2 + s2.^2; % a column vector
    
    % solve the linear system of equations
    cu = b1./a; %a.\b1;
    cv = b2./a; %a.\b2;
    
    % devectorize
    CU = reshape(cu, r, r);
    CV = reshape(cv, r, r);
    
    % Do the similarity transforms
    BU = U*CU*U';
    BV = V*CV*V';
    
    % %% debug
    %
    % norm(R*R'*BU + BU*R*R' - E, 'fro');
    % norm((Sigma.^2)*CU + CU*(Sigma.^2) - E_mod, 'fro');
    % norm(a.*cu - b1, 'fro');
    %
    % norm(R'*R*BV + BV*R'*R - F, 'fro');
    %
    % BU1 = lyap(R*R', - E);
    % norm(R*R'*BU1 + BU1*R*R' - E, 'fro');
    %
    % BV1 = lyap(R'*R, - F);
    % norm(R'*R*BV1 + BV1*R'*R - F, 'fro');
    %
    % % as accurate as the lyap
    % norm(BU - BU1, 'fro')
    % norm(BV - BV1, 'fro')
end



function[Omega1, Omega2] = coupled_lyap(R, E, F)
    % We intent to solve the coupled system of Lyapunov equations
    %
    % RR^T Omega1 + Omega1 RR^T  - R Omega2 R^T = E
    % R^T R Omega2 + Omega1 R^T R  - R^T Omega2 R = F
    %
    % Below is an efficient implementation
    
    [U, Sigma, V] = svd(R);
    E_mod = U'*E*U;
    F_mod = V'*F*V;
    b1 = E_mod(:);
    b2 = F_mod(:);
    
    r = size(Sigma, 1);
    sig = diag(Sigma); % all the singular values in a vector
    sig1 = sig*ones(1, r); % columns repeat
    sig1t = sig1'; % rows repeat
    s1 = sig1(:);
    s2 = sig1t(:);
    
    % The block elements
    a =  s1.^2 + s2.^2; % a column vector
    c = s1.*s2;
    
    % Solve directly using the formula
    % A = diag(a);
    % C = diag(c);
    % Y1_sol = (A*(C\A) - C) \ (b2 + A*(C\b1));
    % Y2_sol = A\(b2 + C*Y1_sol);
    
    Y1_sol = (b2 + (a./c).*b1) ./ ((a.^2)./c - c);
    Y2_sol = (b2 + c.*Y1_sol)./a;
    
    % devectorize
    Omega1 = reshape(Y1_sol, r, r);
    Omega2 = reshape(Y2_sol, r, r);
    
    % Do the similarity transforms
    Omega1 = U*Omega1*U';
    Omega2 = V*Omega2*V';
    
    % %% debug whether we have the right solution
    % norm(R*R'*Omega1 + Omega1*R*R'  - R*Omega2*R' - E, 'fro')
    % norm(R'*R*Omega2 + Omega2*R'*R  - R'*Omega1*R - F, 'fro')
end