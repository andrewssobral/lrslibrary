function [U,output] = cpd_nls(T,U0,options)
%CPD_NLS CPD by nonlinear least squares.
%   [U,output] = cpd_nls(T,U0) computes the factor matrices U{1}, ..., U{N}
%   belonging to a canonical polyadic decomposition of the N-th order
%   tensor T by minimizing 0.5*frob(T-cpdgen(U))^2. The algorithm is
%   initialized with the factor matrices U0{n}. The structure output
%   returns additional information:
%
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   cpd_nls(T,U0,options) may be used to set the following options:
%
%      options.Algorithm =   - The desired optimization method.
%      [@nls_gncgs| ...
%       {@nls_gndl}|@nls_lm]
%      options.LargeScale    - If true, the Gauss-Newton or Levenberg-
%      = sum(size(T))*R>1e2    Marquardt steps are computed using a
%                              preconditioned conjugate gradient algorithm.
%                              Otherwise, a direct solver is used.
%      options.M =           - The preconditioner to use when
%      [{'block-Jacobi'}|...   options.LargeScale is true.
%       'block-SSOR'|false]
%      options.PlaneSearch = - A function handle to the desired CPD plane
%      [{false},@cpd_eps]      search algorithm. Only applicable to
%                              algorithms with trust region globalization.
%      options.<...>         - Parameters passed to the selected method,
%                              e.g., options.TolFun, options.TolX and
%                              options.PlaneSearchOptions. See also help
%                              [options.Algorithm].
%
%   See also cpd_minf.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Optimization-based
%       algorithms for tensor decompositions: canonical polyadic
%       decomposition, decomposition in rank-(Lr,Lr,1) terms and a new
%       generalization," SIAM J. Opt., 2013.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables," SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.

% Format the tensor T.
T = fmt(T,true);
if isstruct(T), size_tens = T.size;
else size_tens = [size(T) ones(1,length(U0)-ndims(T))]; end

% Check the initial factor matrices U0.
N = length(U0);
R = size(U0{1},2);
if any(cellfun('size',U0,2) ~= R)
    error('cpd_nls:U0','size(U0{n},2) should be the same for all n.');
end

% Check the options structure.
isfunc = @(f)isa(f,'function_handle');
xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
if nargin < 3, options = struct; end
if ~isfield(options,'Algorithm')
    funcs = {@nls_gndl,@nls_gncgs,@nls_lm};
    options.Algorithm = funcs{find(cellfun(xsfunc,funcs),1)};
end
if ~isfield(options,'CGMaxIter'), options.CGMaxIter = 10; end
if ~isfield(options,'Display'), options.Display = 0; end
if ~isfield(options,'JHasFullRank'), options.JHasFullRank = false; end
if ~isfield(options,'LargeScale')
    options.LargeScale = sum(size_tens)*R > 1e2;
end
if ~isfield(options,'M'), options.M = 'block-Jacobi'; end
if ~isfield(options,'TolLargeScale'), options.TolLargeScale = 0.02; end

% Adapt line/plane search if it is a CPD line/plane search.
if isfield(options,'LineSearch') && isfunc(options.LineSearch) && ...
   ~isempty(strfind(func2str(options.LineSearch),'cpd_'))
    linesearch = options.LineSearch;
    options.LineSearch = @ls;
end
if isfield(options,'PlaneSearch') && isfunc(options.PlaneSearch) && ...
   ~isempty(strfind(func2str(options.PlaneSearch),'cpd_'))
    planesearch = options.PlaneSearch;
    options.PlaneSearch = @ps;
end

% Call the optimization method.
cache.offset = cumsum([1 cellfun(@numel,U0(:).')]);
cache.T2 = frob(T,'squared');
if options.LargeScale, dF.JHJx = @JHJx; else dF.JHJ = @JHJ; end
dF.JHF = @grad;
switch options.M
    case 'block-SSOR', dF.M = @M_blockSSOR;
    case 'block-Jacobi', dF.M = @M_blockJacobi;
    case 'Jacobi', dF.M = @M_Jacobi;
    otherwise, if isa(options.M,'function_handle'), dF.M = options.M; end
end
state(U0,true);
[U,output] = options.Algorithm(@objfun,dF,U0(:).',options);
output.Name = func2str(options.Algorithm);

function state(z,firstrun)

    if nargin == 2 && firstrun
        % Store the fraction of known elements.
        if isstruct(T) && T.incomplete
            cache.scale = length(T.val)./prod(T.size);
        end
    end

    % Cache the factor matrices' Gramians.
    cache.UHU = zeros(N,R*R);
    for n = 1:N
        tmp = conj(z{n}'*z{n});
        cache.UHU(n,:) = tmp(:);
    end
    
    % Optionally cache the inverses of the Gramians for the preconditioner.
    % In a faster language, this should be the Cholesky factor instead.
    if ischar(options.M) || isa(options.M,'function_handle')
        cache.invW = cell(1,N);
        for n = 1:N
            tmp = cache.UHU([1:n-1 n+1:N],:);
            if N > 2, tmp = prod(tmp,1); end
            cache.invW{n} = inv(reshape(tmp,[R R]));
        end
    end
    
end

function fval = objfun(z)
    
    % CPD objective function.
    isincomplete = isstruct(T) && T.incomplete;
    issparse = isstruct(T) && T.sparse;
    if ~isincomplete || length(T.ind)/prod(T.size) > options.TolLargeScale
        fval = z{1}*kr(z(end:-1:2)).';
        if isincomplete, fval = fval(T.ind)-T.val;
        elseif issparse
            if ~isempty(T.ind), fval(T.ind) = fval(T.ind)-T.val; end
        else fval = fval-reshape(T,size(fval));
        end
    else
        fval = -T.val;
        for r = 1:R
            tmp = z{1}(T.sub{1},r);
            for n = 2:length(z), tmp = tmp.*z{n}(T.sub{n},r); end
            fval = fval+tmp;
        end
    end
    if isincomplete
        E = T; E.val = fval;
        if ~isempty(E.matrix)
            E.matrix = sparse(double(E.sub{1}), ...
                double(1+idivide(E.ind-1,int64(size(E.matrix,1)))), ...
                double(fval),size(E.matrix,1),size(E.matrix,2));
        end
        cache.residual = E;
    else
        cache.residual = reshape(fval,size_tens);
    end
    fval = 0.5*(fval(:)'*fval(:));
    
end

function grad = grad(z)
    
    % CPD scaled conjugate cogradient.
    state(z);
    E = cache.residual;
    offset = cache.offset;
    grad = nan(offset(end)-1,1);
    for n = 1:length(z)
        tmp = full(mtkrprod(E,z,n));
        grad(offset(n):offset(n+1)-1) = tmp(:);
    end
    
end

function JHJ = JHJ(z)
    
    % CPD Jacobian's Gramian.
    UHU = cache.UHU;
    JHJ = zeros(cache.offset(end)-1);
    for n = 1:N
        idxn = cache.offset(n):cache.offset(n+1)-1;
        Wn = reshape(prod(UHU([1:n-1 n+1:N],:),1),[R R]);
        JHJ(idxn,idxn) = kron(Wn,eye(size_tens(n)));
        for m = n+1:N
            idxm = cache.offset(m):cache.offset(m+1)-1;
            Wnm = reshape(prod(UHU([1:n-1 n+1:m-1 m+1:N],:),1),[R R]);
            JHJnm = bsxfun(@times,reshape(z{n},[size_tens(n) 1 1 R]), ...
                    reshape(conj(z{m}),[1 size_tens(m) R 1]));
            JHJnm = bsxfun(@times,JHJnm,reshape(Wnm,[1 1 R R]));
            JHJnm = permute(JHJnm,[1 3 2 4]);
            JHJnm = reshape(JHJnm,[size_tens(n)*R size_tens(m)*R]);
            JHJ(idxn,idxm) = JHJnm;
            JHJ(idxm,idxn) = JHJnm';
        end
    end
    
    % If incomplete, approximate the effect of missing entries.
    if isstruct(T) && T.incomplete
        JHJ = JHJ*cache.scale;
    end
    
end

function y = JHJx(z,x)
    
    % CPD fast Jacobian's Gramian vector product.
    % Ignores the fact that the tensor might be incomplete.
    offset = cache.offset;
    UHU = cache.UHU;
    XHU = zeros(R,R,N);
    y = nan(offset(end)-1,1);
    for n = 1:N
        Wn = UHU([1:n-1 n+1:N],:);
        if N > 2, Wn = prod(Wn,1); end
        tmp = reshape(x(offset(n):offset(n+1)-1),[],R);
        XHU(:,:,n) = conj(tmp'*z{n});
        y(offset(n):offset(n+1)-1) = tmp*reshape(Wn,[R R]);
    end
    for n = 1:N-1
        idxn = offset(n):offset(n+1)-1;
        Wn = zeros(R);
        for m = n+1:N
            idxm = offset(m):offset(m+1)-1;
            if N == 2
                Wn = Wn+XHU(:,:,m);
                JHJmnx = z{m}*XHU(:,:,n);
            else
                Wnm = UHU([1:n-1 n+1:m-1 m+1:N],:);
                if N > 3, Wnm = prod(Wnm,1); end
                Wnm = reshape(Wnm,[R R]);
                Wn = Wn+Wnm.*XHU(:,:,m);
                JHJmnx = z{m}*(Wnm.*XHU(:,:,n));
            end
            y(idxm) = y(idxm)+JHJmnx(:);
        end
        JHJnx = z{n}*Wn;
        y(idxn) = y(idxn)+JHJnx(:);
    end
    
    % If incomplete, approximate the effect of missing entries.
    if isstruct(T) && T.incomplete
        y = y*cache.scale;
    end
    
end

function x = M_blockJacobi(~,b)

    % Solve M*x = b, where M is a block-diagonal approximation for JHJ.
    % Equivalent to simultaneous ALS updates for each of the factors.
    x = nan(size(b));
    for n = 1:length(cache.offset)-1
        idx = cache.offset(n):cache.offset(n+1)-1;
        tmp = reshape(b(idx),[],size(cache.invW{1},1))*cache.invW{n};
        x(idx) = tmp(:);
    end
    
    % If incomplete, approximate the effect of missing entries.
    if isstruct(T) && T.incomplete
        x = x/cache.scale;
    end
    
end

function x = M_blockSSOR(z,b)
    
    % Solve Mx = b, where M is a block-Symmetric Successive Overrelaxation
    % preconditioner.
    % x = inv(D)*(U+D)*b
    B = cell(size(z));
    UHU = cache.UHU;
    BHU = nan(size(UHU));
    for n = 1:N
        B{n} = b(cache.offset(n):cache.offset(n+1)-1);
        BHU(:,:,n) = B{n}'*z{n};
    end
    X = B;
    for n = 1:N-1
        Wsum = zeros(R);
        for m = n+1:N
            Wnm = conj(reshape(prod(UHU([1:n-1 n+1:m-1 m+1:N],:),1),R,[]));
            Wsum = Wsum+Wnm.*conj(BHU(:,:,m));
        end
        Wn = conj(reshape(prod(UHU([1:n-1 n+1:N],:),1),[R R]));
        X{n} = X{n}+(z{n}*Wsum)/Wn;
    end
    % x = (L+D)*x
    B = X;
    for n = 1:N
        Wn = conj(reshape(prod(UHU([1:n-1 n+1:N],:),1),[R R]));
        BHU(:,:,n) = B{n}'*z{n};
        X{n} = B{n}*Wn;
    end
    for n = 2:N
        Wsum = zeros(R);
        for m = 1:n-1
            Wnm = conj(reshape(prod(UHU([1:n-1 n+1:m-1 m+1:N],:),1),R,[]));
            Wsum = Wsum+Wnm.*conj(BHU(:,:,m));
        end
        X{n} = X{n}+z{n}*Wsum;
    end
    x = cell2mat(cellfun(@(x)x(:),X(:),'UniformOutput',false));
    
    % If incomplete, approximate the effect of missing entries.
    if isstruct(T) && T.incomplete
        x = x/cache.scale;
    end
    
end

function [alpha,output] = ls(~,~,z,p,state,options)
    state.UHU = cache.UHU;
    state.T2 = cache.T2;
    [alpha,output] = linesearch(T,z,p,state,options);
end

function [alpha,output] = ps(~,~,z,p,q,state,options)
    state.UHU = cache.UHU;
    state.T2 = cache.T2;
    [alpha,output] = planesearch(T,z,p,q,state,options);
end

end
