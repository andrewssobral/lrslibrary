function [U,output] = btd_minf(T,U0,options)
%BTD_MINF BTD by unconstrained nonlinear optimization.
%   [U,output] = btd_minf(T,U0) computes R terms U{r} corresponding to a
%   block term decomposition of the N-th order tensor T by minimizing
%   0.5*frob(T-btdgen(U))^2. Each term U{r} is a cell array of N factor
%   matrices U{r}{n}, followed by a core tensor U{r}{N+1}. The algorithm is
%   initialized with the terms matrices U0{r}. The structure output returns
%   additional information:
%
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   btd_minf(T,U0,options) may be used to set the following options:
%
%      options.Algorithm =     - The desired optimization method.
%      [{@minf_lbfgsdl}|...
%       @minf_lbfgs|@minf_ncg]
%      options.<...>           - Parameters passed to the selected method,
%                                e.g., options.TolFun, options.TolX and
%                                options.LineSearchOptions. See also help
%                                [options.Algorithm].
%
%   See also btd_nls.

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

% Check the options structure.
isfunc = @(f)isa(f,'function_handle');
xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
if nargin < 3, options = struct; end
if ~isfield(options,'Algorithm')
    funcs = {@minf_lbfgsdl,@minf_lbfgs,@minf_ncg};
    options.Algorithm = funcs{find(cellfun(xsfunc,funcs),1)};
end
if ~isfield(options,'Display'), options.Display = 0; end
if ~isfield(options,'TolFun'), options.TolFun = 1e-12; end
if ~isfield(options,'TolLargeScale'), options.TolLargeScale = 0.02; end

% Call the optimization method.
cache.offset = cell2mat(cellfun(@(t)cellfun(@numel,t(:).'),U0(:).', ...
    'UniformOutput',false));
cache.offset = cumsum([1 cache.offset]);
cache.T2 = frob(T,'squared');
[U,output] = options.Algorithm(@objfun,@grad,U0(:).',options);
output.Name = func2str(options.Algorithm);

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
    N = length(z{1})-1;
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

end
