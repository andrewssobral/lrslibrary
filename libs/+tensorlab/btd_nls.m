function [U,output] = btd_nls(T,U0,options)
%BTD_NLS BTD by nonlinear least squares.
%   [U,output] = btd_nls(T,U0) computes R terms U{r} corresponding to a
%   block term decomposition of the N-th order tensor T by minimizing
%   0.5*frob(T-btdgen(U))^2. Each term U{r} is a cell array of N factor
%   matrices U{r}{n}, followed by a core tensor U{r}{N+1}. The algorithm is
%   initialized with the terms matrices U0{r}. The structure output returns
%   additional information:
%
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   btd_nls(T,U0,options) may be used to set the following options:
%
%      options.Algorithm =   - The desired optimization method.
%      [@nls_gncgs| ...
%       {@nls_gndl}|@nls_lm]
%      options.M =           - The preconditioner to use when
%      [{'block-Jacobi'}|...   options.LargeScale is true.
%       false]
%      options.<...>         - Parameters passed to the selected method,
%                              e.g., options.TolFun, options.TolX and
%                              options.PlaneSearchOptions. See also help
%                              [options.Algorithm].
%
%   See also btd_minf.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Unconstrained
%       optimization of real functions in complex variables," SIAM J. Opt.,
%       Vol. 22, No. 3, 2012, pp. 879-898.

% Format the tensor T.
T = fmt(T,true);
if isstruct(T), size_tens = T.size;
else size_tens = [size(T) ones(1,length(U0)-ndims(T))]; end
N = length(U0{1})-1;
R = length(U0);

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
if ~isfield(options,'M'), options.M = 'block-Jacobi'; end
if ~isfield(options,'TolLargeScale'), options.TolLargeScale = 0.02; end

% Call the optimization method.
cache.offset = cell2mat(cellfun(@(t)cellfun(@numel,t(:).'),U0(:).', ...
    'UniformOutput',false));
cache.offset = cumsum([1 cache.offset]);
cache.T2 = frob(T,'squared');
dF.JHJx = @JHJx;
dF.JHF = @grad;
switch options.M
    case 'block-Jacobi', dF.M = @M_blockJacobi;
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
    [idx,jdx,kdx] = ndgrid(1:R,1:N,1:R);
    cache.UHU = ...
        arrayfun(@(i,n,j)z{i}{n}'*z{j}{n}, ...
        idx,jdx,kdx,'UniformOutput',false);
    
    % Optionally cache some results for the block-Jacobi preconditioner.
    if ischar(options.M) || isa(options.M,'function_handle')
        [idx,jdx] = ndgrid(1:R,1:N);
        UHU = cache.UHU;
        cache.invSKS = arrayfun( ...
            @(r,n)inv(mtkronprod(z{r}{end},UHU(r,:,r),n)* ...
            conj(reshape(permute(z{r}{end},[1:n-1 n+1:N n]), ...
            [],size(z{r}{end},n)))),idx,jdx,'UniformOutput',false);
        cache.invUHU = arrayfun( ...
            @(r,n)inv(UHU{r,n,r}),idx,jdx,'UniformOutput',false);
    end

end

function fval = objfun(z)
    
    % BTD objective function.
    isincomplete = isstruct(T) && T.incomplete;
    issparse = isstruct(T) && T.sparse;
    if ~isincomplete || length(T.ind)/prod(T.size) > options.TolLargeScale
        fval = z{1}{1}*mtkronprod(z{1}{end},z{1}(1:end-1),1,'H');
        for r = 2:length(z)
            fval = fval+z{r}{1}*mtkronprod(z{r}{end},z{r}(1:end-1),1,'H');
        end
        if isincomplete, fval = fval(T.ind)-T.val;
        elseif issparse
            if ~isempty(T.ind), fval(T.ind) = fval(T.ind)-T.val; end
        else fval = fval-reshape(T,size(fval));
        end
    else
        fval = -T.val;
        for r = 1:length(z)
            size_core = cellfun('size',z{r}(1:end-1),2);
            idx = cell(1,length(size_core));
            S = z{r}{end};
            for i = 1:numel(S)
                [idx{:}] = ind2sub(size_core,i);
                tmp = S(idx{:})*z{r}{1}(T.sub{1},idx{1});
                for n = 2:length(size_core)
                    tmp = tmp.*z{r}{n}(T.sub{n},idx{n});
                end
                fval = fval+tmp;
            end
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
    
    % BTD scaled conjugate cogradient.
    state(z);
    E = cache.residual;
    offset = cache.offset;
    grad = nan(offset(end)-1,1);
    cnt = 1;
    for r = 1:length(z)
        V = z{r}(1:N);
        S = conj(z{r}{end});
        for n = 1:N
            tmp = full(mtkronprod(E,V,n))* ...
                reshape(permute(S,[1:n-1 n+1:N n]),[],size(S,n));
            grad(offset(cnt):offset(cnt+1)-1) = tmp(:);
            cnt = cnt+1;
        end
        tmp = full(mtkronprod(E,V,0));
        grad(offset(cnt):offset(cnt+1)-1) = tmp;
        cnt = cnt+1;
    end
    
end

function y = JHJx(z,x)
    
    % BTD fast Jacobian's Gramian vector product.
    % Ignores the fact that the tensor might be incomplete.
    offset = cache.offset;
    UHU = cache.UHU;
    [idx,jdx,kdx] = ndgrid(1:R,1:N,1:R);
    x = deserialize(x);
    XHU = arrayfun(@(i,n,j)x{i}{n}'*z{j}{n}, ...
        idx,jdx,kdx,'UniformOutput',false);
    y = nan(offset(end)-1,1);
    cnt = 1;
    for ri = 1:R
        
        % Factor matrices.
        for ni = 1:N
            idx = offset(cnt):offset(cnt+1)-1;
            Sri = permute(z{ri}{end},[1:ni-1 ni+1:N ni]);
            Sri = conj(reshape(Sri,[],size(Sri,N)));
            for rj = 1:R
                Srj = z{rj}{end};
                tmp = mtkronprod(x{rj}{end},UHU(rj,:,ri),ni);
                for nj = [1:ni-1 ni+1:N]
                    proj = UHU(rj,:,ri);
                    proj{nj} = XHU{rj,nj,ri};
                    tmp = tmp+mtkronprod(Srj,proj,ni);
                end
                tmp = z{rj}{ni}*(tmp*Sri);
                tmp = tmp+x{rj}{ni}* ...
                    (mtkronprod(z{rj}{end},UHU(rj,:,ri),ni)*Sri);
                if rj == 1, y(idx) = tmp(:);
                else y(idx) = y(idx)+tmp(:); end
            end
            cnt = cnt+1;
        end
        
        % Core tensor.
        idx = offset(cnt):offset(cnt+1)-1;
        for rj = 1:R
            Srj = z{rj}{end};
            tmp = mtkronprod(x{rj}{end},UHU(rj,:,ri),0);
            for nj = 1:N
                proj = UHU(rj,:,ri);
                proj{nj} = XHU{rj,nj,ri};
                tmp = tmp+mtkronprod(Srj,proj,0);
            end
            if rj == 1, y(idx) = tmp(:);
            else y(idx) = y(idx)+tmp(:); end
        end
        cnt = cnt+1;
        
    end
    
    % If incomplete, approximate the effect of missing entries.
    if isstruct(T) && T.incomplete
        y = y*cache.scale;
    end
    
end

function x = M_blockJacobi(~,b)

    % Solve M*x = b, where M is a block-diagonal approximation for JHJ.
    x = nan(size(b));
    offset = cache.offset;
    invSKS = cache.invSKS;
    invUHU = cache.invUHU;
    cnt = 1;
    for r = 1:R
        for n = 1:N
            idx = offset(cnt):offset(cnt+1)-1;
            tmp = reshape(b(idx),[],size(invSKS{r,n},1))*invSKS{r,n};
            x(idx) = tmp(:);
            cnt = cnt+1;
        end
        idx = offset(cnt):offset(cnt+1)-1;
        size_core = cellfun('size',invUHU(r,:),1);
        x(idx) = mtkronprod(reshape(b(idx),size_core),invUHU(r,:),0);
        cnt = cnt+1;
    end
    
    % If incomplete, approximate the effect of missing entries.
    if isstruct(T) && T.incomplete
        x = x/cache.scale;
    end
    
end

function x = deserialize(x)
    off = cache.offset; tmp = x; cnt = 1;
    x = cell(1,R);
    for r = 1:R
        x{r} = cell(1,N+1);
        for n = 1:N
            x{r}{n} = reshape(tmp(off(cnt):off(cnt+1)-1),size(U0{r}{n}));
            cnt = cnt+1;
        end
        x{r}{end} = reshape(tmp(off(cnt):off(cnt+1)-1),size(U0{r}{end}));
        cnt = cnt+1;
    end
end

end
