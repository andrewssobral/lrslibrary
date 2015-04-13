function [U,output] = cpd_minf(T,U0,options)
%CPD_MINF CPD by unconstrained nonlinear optimization.
%   [U,output] = cpd_minf(T,U0) computes the factor matrices U{1}, ...,
%   U{N} belonging to a canonical polyadic decomposition of the N-th order
%   tensor T by minimizing 0.5*frob(T-cpdgen(U))^2. The algorithm is
%   initialized with the factor matrices U0{n}. The structure output
%   returns additional information:
%
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   cpd_minf(T,U0,options) may be used to set the following options:
%
%      options.Algorithm =     - The desired optimization method.
%      [{@minf_lbfgsdl}|...
%       @minf_lbfgs|@minf_ncg]
%      options.LineSearch =    - A function handle to the desired line
%      [{'auto'},@ls_mt|...      search algorithm. If the line search
%       @cpd_aels|@cpd_els|...   method has the prefix 'cpd_', it is
%       @cpd_lsb]                modified to be compatible with the
%                                optimization algorithm. Only applicable to
%                                algorithms with line search globalization.
%      options.PlaneSearch =   - A function handle to the desired CPD plane
%      [{false},@cpd_eps]        search algorithm. Only applicable to
%                                algorithms with trust region
%                                globalization.
%      options.<...>           - Parameters passed to the selected method,
%                                e.g., options.TolFun, options.TolX and
%                                options.LineSearchOptions. See also help
%                                [options.Algorithm].
%
%   See also cpd_nls.

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
R = size(U0{1},2);
if any(cellfun('size',U0,2) ~= R)
    error('cpd_minf:U0','size(U0{n},2) should be the same for all n.');
end

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

% Adapt line/plane search if it is a CPD line/plane search.
if isfield(options,'LineSearch') && ...
   ~isempty(strfind(func2str(options.LineSearch),'cpd_'))
    linesearch = options.LineSearch;
    options.LineSearch = @ls;
end
if isfield(options,'PlaneSearch') && ...
   ~isempty(strfind(func2str(options.PlaneSearch),'cpd_'))
    planesearch = options.PlaneSearch;
    options.PlaneSearch = @ps;
end

% Call the optimization method.
cache.offset = cumsum([1 cellfun(@numel,U0(:).')]);
cache.T2 = frob(T,'squared');
[U,output] = options.Algorithm(@objfun,@grad,U0(:).',options);
output.Name = func2str(options.Algorithm);

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
    E = cache.residual;
    offset = cache.offset;
    grad = nan(offset(end)-1,1);
    for n = 1:length(z)
        tmp = full(mtkrprod(E,z,n));
        grad(offset(n):offset(n+1)-1) = tmp(:);
    end
    
end

function [alpha,output] = ls(~,~,z,p,state,options)
    state.T2 = cache.T2;
    [alpha,output] = linesearch(T,z,p,state,options);
end

function [alpha,output] = ps(~,~,z,p,q,state,options)
    state.T2 = cache.T2;
    [alpha,output] = planesearch(T,z,p,q,state,options);
end

end
