function [sol,output] = sdf_minf(model,options)
%SDF_MINF Structured data fusion by unconstrained nonlinear optimization.
%   [sol,output] = sdf_minf(model) solves the data fusion problem described
%   by the structure model and returns the solution as the structure sol.
%   The model describes a data fusion problem with three fields:
%
%      model.variables
%
%         A structure or cell array of initializations for the variables in
%         the data fusion problem. Each field or cell in model.variables 
%         represents a variable. A variable may be an array (such as
%         scalars, vectors, matrices and tensors) or a (nested) cell array
%         of arrays.
%
%      model.factors
%
%         A structure or cell array of factors. Each field or cell in
%         model.factors represents a factor. A factor is described by a
%         cell array of subfactors. After each subfactor has been
%         generated, the factor is constructed as the cell2mat of its
%         subfactors. A subfactor is a cell array in which the first
%         element is a reference to a variable, and the following elements
%         represent a sequence of transformations of that variable. If
%         model.variables is a structure, then a reference to a variable is
%         a string corresponding to the field name of the desired variable.
%         If model.variables is a cell array, then a reference to a
%         variable is the index of the cell containing the desired
%         variable. Transformations are supplied by functions which
%         describe their (linearized) behaviour. For example, the
%         transformation @struct_inv computes the factor as the matrix
%         inverse of the selected variable. All transformations must have
%         the function signature struct_mytrans(z,task). Instead of
%         supplying a reference to a variable, it is also possible to
%         supply a constant factor in the form of an array.
%         
%      model.factorizations
%         
%         A structure of data sets to jointly factorize. Each field in
%         model.factorizations represents a factorization. A factorization
%         is described by a structure containing two fields. The first
%         field has the field name 'data' and contains the (dense, sparse
%         or incomplete) array which is to be factorized. If the array
%         contains many zeros or a NaN, it is internally converted to a
%         sparse or incomplete tensor with fmt. The second field's 
%         field name designates the type of factorization to compute of the
%         data set. Currently, two types are supported:
%
%         1. Canonical polyadic decomposition: a field 'cpd' should contain
%            a cell array of references to factors in model.factors. The
%            nth reference in the cell array corresponds to the nth factor
%            matrix of the CPD.
%
%         2. Block term decomposition: a field 'btd' should contain a cell
%            array of terms. Each term is itself a cell array containing
%            references to factors in model.factors. The nth reference in a
%            term corresponds to the nth factor matrix of that term. The
%            (N+1)th reference, where N is the number of dimensions of the
%            data set, in a term corresponds to the core tensor of that
%            term.
%
%         Additionally, two types of regularization are available:
%
%         3. L2 regularization: a field 'regL2' should contain a cell array
%            of references to factors in model.factors. This appends a term
%            to the objective function of the form 0.5*norm(F(:)-D(:),2)^2,
%            where F(:) and D(:) are the serialized factors in regL2 and
%            the data field respectively. The data field may be omitted, in
%            which case it is an all-zero vector.
%
%         4. L1 regularization: a field 'regL1' should contain a cell array
%            of references to factors in model.factors. This appends a term
%            to the objective function which is a smooth approximation of
%            0.5*norm(F(:)-D(:),1), where F(:) and D(:) are the serialized
%            factors in regL1 and the data field respectively. The data
%            field may be omitted, in which case it is an all-zero vector.
%
%   The algorithm proceeds to compute the joint factorization of the data
%   sets provided by models.factorizations by minimizing the sum of square
%   magnitude residuals of these factorizations. The output sol contains
%   the optimized variables in sol.variables and the corresponding factors
%   in sol.factors.
%
%   The structure output returns additional information:
%
%      output.Name  - The name of the selected algorithm.
%      output.<...> - The output of the selected algorithm.
%
%   sdf_minf(model,options) may be used to set the following options:
%
%      options.Algorithm =     - The desired optimization method.
%      [{@minf_lbfgsdl}|...
%       @minf_lbfgs|@minf_ncg]
%      options.RelWeights =    - By supplying relative weights, the weights
%      ones(1,F)                 options.Weights are computed as follows:
%                                options.Weights(f) = options.RelWeights(f)
%                                /(sum(options.RelWeights)*numel(data_f))
%                                for each of the F factorizations in the
%                                model.
%      options.Weights         - The weight of the F factorizations in the
%                                data fusion model. The SDF objective
%                                function is \sum_f 0.5*options.Weights(f)*
%                                frob(model_f-data_f)^2, where model_f and
%                                data_f are the fth factorization and 
%                                corresponding tensor. By default, weights
%                                are provided by options.RelWeights, but
%                                options.Weights has precedence if
%                                supplied.
%      options.<...>           - Parameters passed to the selected method,
%                                e.g., options.TolFun, options.TolX.
%                                See also help [options.Algorithm].
%
%   See also sdf_nls.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

% Check the options structure.
if nargin < 2, options = struct; end
isfunc = @(f)isa(f,'function_handle');
xsfunc = @(f)isfunc(f)&&exist(func2str(f),'file');
if ~isfield(options,'Algorithm')
    funcs = {@minf_lbfgsdl,@minf_lbfgs,@minf_ncg};
    options.Algorithm = funcs{find(cellfun(xsfunc,funcs),1)};
end
if ~isfield(options,'MaxIter'), options.MaxIter = 5000; end
if ~isfield(options,'Display'), options.Display = 0; end
if ~isfield(options,'TolFun'), options.TolFun = 1e-12; end
if ~isfield(options,'TolLargeScale'), options.TolLargeScale = 0.02; end

% Convert model to internal format.
if ~iscell(model.variables)
    cache.variables.names = fieldnames(model.variables);
    model.variables = reshape(struct2cell(model.variables),1,[]);
else
    model.variables = model.variables(:).';
end
if ~iscell(model.factors)
    cache.factors.names = fieldnames(model.factors);
    model.factors = reshape(struct2cell(model.factors),1,[]);
else
    model.factors = model.factors(:).';
end
if ~iscell(model.factorizations)
    cache.factorizations.names = fieldnames(model.factorizations);
    model.factorizations = reshape(struct2cell(model.factorizations),1,[]);
else
    model.factorizations = model.factorizations(:).';
end
for I = 1:length(model.factors)
    if ~iscell(model.factors{I})
        model.factors{I} = model.factors(I);
    end
    if any(cellfun(@(f)isa(f,'function_handle'),model.factors{I}))
        model.factors{I} = model.factors(I);
    else
        for J = find(~cellfun(@iscell,model.factors{I}(:).'))
            model.factors{I}{J} = model.factors{I}(J);
        end
    end
end

% Set functions for saving state, computing the objective function,
% gradient and fast matrix-vector products. New models need only implement
% these three functions.
for I = 1:length(model.factorizations)
    fn = fieldnames(model.factorizations{I});
    model.factorizations{I}.type = find(~strcmp('data',fn));
    model.factorizations{I}.type = fn{model.factorizations{I}.type};
    model.factorizations{I}.factors = ...
        model.factorizations{I}.(model.factorizations{I}.type);
    model.factorizations{I} = ...
        rmfield(model.factorizations{I},model.factorizations{I}.type);
    switch model.factorizations{I}.type
        case 'cpd'
            model.factorizations{I}.state = @state_cpd;
            model.factorizations{I}.objfun = @objfun_cpd;
            model.factorizations{I}.grad = @grad_cpd;
        case 'btd'
            model.factorizations{I}.state = @state_btd;
            model.factorizations{I}.objfun = @objfun_btd;
            model.factorizations{I}.grad = @grad_btd;
        case 'regL2'
            model.factorizations{I}.state = @state_regL2;
            model.factorizations{I}.objfun = @objfun_regL2;
            model.factorizations{I}.grad = @grad_regL2;
        case 'regL1'
            model.factorizations{I}.state = @state_regL1;
            model.factorizations{I}.objfun = @objfun_regL1;
            model.factorizations{I}.grad = @grad_regL1;
        otherwise
            error('sdf_minf:model','Model %s is not supported', ...
                model.factorizations{I}.type);
    end
end

% Fill in constants and dereference pointers from factors to variables.
cache.factors.isconst = cell(1,length(model.factors));
for I = 1:length(model.factors)
    cache.factors.isconst{I} = ...
        cellfun(@(f)isnumeric(f{1})&&~isscalar(f{1}),model.factors{I});
    for J = 1:numel(model.factors{I})
        if cache.factors.isconst{I}(J)
            const = model.factors{I}{J}{1};
            for K = 2:length(model.factors{I}{J})
                const = model.factors{I}{J}{K}(const);
            end
            model.factors{I}{J} = {const};
        elseif ischar(model.factors{I}{J}{1})
            model.factors{I}{J}{1} = ...
                find(strcmp(model.factors{I}{J}{1},cache.variables.names));
        end
    end
end

% Dereference and flatten pointers from factorizations to factors.
for I = 1:length(model.factorizations)
	model.factorizations{I}.factors = ...
        deref(model.factorizations{I}.factors);
end
function array = deref(array)
    for i = 1:length(array)
        if iscell(array{i})
            array{i} = deref(array{i});
        elseif ischar(array{i})
            array{i} = find(strcmp(array{i},cache.factors.names));
        end
    end
end

% Convert full data to incomplete/sparse format where appropriate.
cache.factorizations.isincomplete = false(1,length(model.factorizations));
cache.factorizations.issparse = false(1,length(model.factorizations));
for I = 1:length(model.factorizations)
    if ~isfield(model.factorizations{I},'data') || ...
       iscell(model.factorizations{I}.data)
        continue;
    end
    model.factorizations{I}.data = ...
        fmt(model.factorizations{I}.data,true);
    if isstruct(model.factorizations{I}.data)
        cache.factorizations.isincomplete(I) = ...
            model.factorizations{I}.data.incomplete;
        cache.factorizations.issparse(I) = ...
            model.factorizations{I}.data.sparse;
    end
end

% Cache expansion of variables into factors.
cache.factors.expanded = cell(size(model.factors));
cache.factors.sequence = cell(size(model.factors));
cache.factors.state = cell(size(model.factors));
for I = 1:length(model.factors)
    cache.factors.sequence{I} = cell(size(model.factors{I}));
    cache.factors.state{I} = cell(size(model.factors{I}));
    for J = 1:numel(model.factors{I})
        cache.factors.sequence{I}{J} = cell(1,length(model.factors{I}{J}));
        cache.factors.state{I}{J} = cell(1,length(model.factors{I}{J})-1);
    end
end
expand(model.variables,'cache');

% Cache factors' structure.
cache.factors.structure = ...
    cellfun(@structure,cache.factors.sequence,'UniformOutput',false);

% Cache offsets for variables and expanded factors.
cache.variables.offset = ...
    cumsum([1 cellfun(@(v)numel(serialize(v)),model.variables)]);
cache.factorizations.serialized = cell(size(model.factorizations));
cache.factorizations.offset = cell(size(model.factorizations));
cache.factorizations.suboffset = cell(size(model.factorizations));
for I = 1:length(model.factorizations)
    cache.factorizations.serialized{I} = ...
        serialize(model.factorizations{I}.factors).';
    cache.factorizations.offset{I} = ...
        cumsum([1 cellfun(@(v)numel(serialize(v)), ...
        cache.factors.expanded(cache.factorizations.serialized{I}))]);
    cum = cumsum([0 cellfun(@(f)sum(serialize( ...
        cellfun(@(s)numel(s{end}),f))),cache.factors.sequence( ...
        cache.factorizations.serialized{I}))]);
    fct = cellfun(@(f)cellfun(@(s)false(size(s{end})),f,'UniformOutput',...
        false),cache.factors.sequence( ...
        cache.factorizations.serialized{I}),'UniformOutput',false);
    cache.factorizations.suboffset{I} = cell(1,length(fct));
    for J = 1:length(fct)
        cache.factorizations.suboffset{I}{J} = cell(size(fct{J}));
        for K = 1:numel(fct{J})
            subfct = fct{J};
            subfct{K} = true(size(subfct{K}));
            subfct = find(cell2mat(subfct));
            if all(subfct(2:end) == subfct(1:end-1)+1)
                subfct = [subfct(1) subfct(end)];
            end
            cache.factorizations.suboffset{I}{J}{K} = subfct+cum(J);
        end
    end
end

% Set (relative) weights.
if isfield(options,'RelWeights') && isfield(options,'Weights')
    warning('sdf_nls:weights',['Both relative and absolute weights ' ...
        'are supplied, proceeding with absolute weights.']);
end
if ~isfield(options,'RelWeights')
    options.RelWeights = ones(1,length(model.factorizations));
end
options.RelWeights = options.RelWeights./sum(options.RelWeights);
if ~isfield(options,'Weights')
    for I = 1:length(model.factorizations)
        if isfield(model.factorizations{I},'data') && ...
           ~iscell(model.factorizations{I}.data)
            if isstruct(model.factorizations{I}.data)
                NUM = numel(model.factorizations{I}.data.val);
            else
                NUM = numel(model.factorizations{I}.data);
            end
        else
            NUM = cache.factorizations.offset{I}(end);
        end
        options.Weights(I) = 2*options.RelWeights(I)/NUM;
    end
end
if length(options.Weights) ~= length(model.factorizations)
    error('sdf_nls:weights',['The number of weights must equal the ' ...
        'number of factorizations.']);
end

% Initialize each factorization's state.
cache.factorizations.state = cell(size(model.factorizations));
for I = 1:length(model.factorizations)
    model.factorizations{I}.state(I,true);
end

% Run optimization algorithm.
[z,output] = options.Algorithm(@objfun,@grad,model.variables,options);
output.Name = func2str(options.Algorithm);

% Return output.
if isfield(cache.variables,'names')
    sol.variables = cell2struct(z(:),cache.variables.names);
else
    sol.variables = z;
end
if isfield(cache.factors,'names')
    sol.factors = cell2struct(reshape(expand(z),[],1),cache.factors.names);
else
    sol.factors = expand(z);
end

function x = expand(z,where)
% Expands the variables into factors stored in the field factors.expanded,
% and stores their computational state in the field factors.state.

    % Save some references for speed.
    factors = model.factors;
    isconst = cache.factors.isconst;
    
    % Expand into cache or into the output variable x.
    if nargin == 2 && ischar(where)
        
        % For each ith factor...
        for i = 1:length(factors)
            % For each jth subfactor...
            for j = 1:numel(factors{i})
                if isconst{i}(j)
                    % Constant subfactor.
                    cache.factors.sequence{i}{j}{1} = factors{i}{j}{1};
                else
                    % For each kth transformation...
                    cache.factors.sequence{i}{j}{1} = z{factors{i}{j}{1}};
                    for k = 2:length(factors{i}{j})
                        [cache.factors.sequence{i}{j}{k}, ...
                            cache.factors.state{i}{j}{k-1}] = ...
                            factors{i}{j}{k}( ...
                            cache.factors.sequence{i}{j}{k-1},[]);
                    end
                end
            end
            if numel(cache.factors.sequence{i}) == 1
                cache.factors.expanded{i} = ...
                    cache.factors.sequence{i}{j}{end};
            else
                cache.factors.expanded{i} = ...
                    cell2mat(cellfun(@(f)f{end}, ...
                    cache.factors.sequence{i},'UniformOutput',false));
            end
        end
        
    else
        
        % For each ith factor...
        x = cell(size(factors));
        for i = 1:length(factors)
            % For each jth subfactor...
            sub = cell(size(factors{i}));
            for j = 1:numel(factors{i})
                if isconst{i}(j)
                    % Constant subfactor.
                    sub{j} = factors{i}{j}{1};
                else
                    % For each kth transformation...
                    sub{j} = z{factors{i}{j}{1}};
                    for k = 2:length(factors{i}{j})
                        sub{j} = factors{i}{j}{k}(sub{j},[]);
                    end
                end
            end
            if numel(sub) == 1, x{i} = sub{1};
            else x{i} = cell2mat(sub);
            end
        end
        
    end
    
end

function y = derivcontract(f,y,l)
% Linearly contract factors using their transformations' Jacobians.

    % Save some references for speed.
    factors = model.factors;
    voffset = cache.variables.offset;
    isconst = cache.factors.isconst;
    sequence = cache.factors.sequence;
    state = cache.factors.state;
    structure = cache.factors.structure;
    serialized = cache.factorizations.serialized{f};
    soffset = cache.factorizations.suboffset{f};

    % Update y with Jacobian-contracted variables.
    for i = 1:length(serialized)
        % For each jth subfactor...
        idx = serialized(i);
        for j = 1:numel(factors{idx})
            
            % Skip this subfactor if it's constant.
            if isconst{idx}(j), continue; end

            % Apply sequence of Jacobian-vector products.
            seq = factors{idx}{j};
            ftr = sequence{idx}{j};
            stt = state{idx}{j};
            if size(soffset{i}{j},2) == 2
                sub = l(soffset{i}{j}(1):soffset{i}{j}(2));
            else
                sub = l(soffset{i}{j});
            end
            if length(seq) > 1
                dim = structure{idx}{j}{end};
                if isnumeric(dim)
                    if ~isempty(dim), sub = reshape(sub,dim); end
                else sub = deserialize(sub,dim);
                end
                for k = length(seq):-1:2
                    task = stt{k-1};
                    task.l = sub;
                    task.r = [];
                    sub = seq{k}(ftr{k-1},task);
                end
                sub = serialize(sub);
            end
            
            % Update y.
            jdx = voffset(seq{1}):voffset(seq{1}+1)-1;
            y(jdx) = y(jdx)+sub;
        
        end
    end
    
end

function z = getfactors(f,x)
% Retrieves factors of the fth factorization, given the factors x.
    
    % Deserialize the factors.
    z = deserialize_local(x,model.factorizations{f}.factors);
    function z = deserialize_local(x,dim)
        z = cell(size(dim));
        for i = 1:numel(dim)
            if iscell(dim{i})
                z{i} = deserialize_local(x,dim{i});
            else
                z{i} = x{dim{i}};
            end
        end
    end
    
end

function fval = objfun(z)
    
    % Expand the variables into factors.
    x = expand(z);
    
    % Compute objective function value.
    fval = 0;
    for i = find(options.Weights(:).' ~= 0)
        fval = fval+options.Weights(i)* ...
            model.factorizations{i}.objfun(i,getfactors(i,x));
    end
    
end

function grad = grad(z)
    
    % Expand the variables into factors.
    expand(z,'cache');
    
    % Let each factorization save intermediate computations in the cache.
    for i = 1:length(model.factorizations)
        model.factorizations{i}.state(i);
    end
    
    % Compute the gradient.
    grad = zeros(cache.variables.offset(end)-1,1);
    for i = find(options.Weights(:).' ~= 0)
        
        % Compute model's gradient.
        tmp = getfactors(i,cache.factors.expanded);
        JHF = model.factorizations{i}.grad(i,tmp);
        
        % Contract the factor matrices into variables.
        if i == 1
            grad = options.Weights(i)*derivcontract(i,grad,JHF);
        else
            grad = grad+options.Weights(i)*derivcontract(i,grad,JHF);
        end
        
    end
    
end

function [z,offset] = deserialize(z,dim,offset)
    if iscell(dim)
        v = z;
        z = cell(size(dim));
        if nargin < 3, offset = 0; end
        for i = 1:numel(z)
            if iscell(dim{i})
                [z{i},offset] = deserialize(v,dim{i},offset);
            else
                n = prod(dim{i}(:));
                z{i} = reshape(v(offset+(1:n)),dim{i});
                offset = offset+n;
            end
        end
    elseif ~isempty(dim)
        z = reshape(z,dim);
    end
end

function z = serialize(z)
    if iscell(z)
        for i = find(cellfun(@iscell,z(:).'))
            z{i} = serialize(z{i});
        end
        s = cellfun(@numel,z(:)); o = [0; cumsum(s)];
        c = z; z = zeros(o(end),1);
        for i = 1:length(s), z(o(i)+(1:s(i))) = c{i}(:); end
    else
        z = z(:);
    end
end

function dim = structure(z)
    if iscell(z)
        dim = cellfun(@size,z,'UniformOutput',false);
        for i = find(cellfun(@iscell,z(:).'))
            dim{i} = structure(z{i});
        end
    else
        dim = size(z);
        if numel(z) == dim(1), dim = []; end
    end
end

% Model: canonical polyadic decomposition ---------------------------------

function state_cpd(~,~)
end

function fval = objfun_cpd(f,z)
    
    % CPD objective function.
    T = model.factorizations{f}.data;
    isincomplete = cache.factorizations.isincomplete(f);
    issparse = cache.factorizations.issparse(f);
    if ~isincomplete || length(T.ind)/prod(T.size) > options.TolLargeScale
        fval = z{1}*kr(z(end:-1:2)).';
        if isincomplete, fval = fval(T.ind)-T.val;
        elseif issparse
            if ~isempty(T.ind), fval(T.ind) = fval(T.ind)-T.val; end
        else fval = fval-reshape(T,size(fval));
        end
    else
        fval = -T.val;
        for r = 1:size(z{1},2)
            tmp = z{1}(T.sub{1},r);
            for n = 2:length(z), tmp = tmp.*z{n}(T.sub{n},r); end
            fval = fval+tmp;
        end
    end
    if isincomplete
        T.val = fval;
        if ~isempty(T.matrix)
            T.matrix = sparse(double(T.sub{1}), ...
                double(1+idivide(T.ind-1,int64(size(T.matrix,1)))), ...
                double(fval),size(T.matrix,1),size(T.matrix,2));
        end
        cache.factorizations.residual{f} = T;
    else
        if issparse, size_tens = T.size;
        else size_tens = size(T); end
        cache.factorizations.residual{f} = reshape(fval,size_tens);
    end
    fval = 0.5*(fval(:)'*fval(:));
    
end

function grad = grad_cpd(f,z)
    
    % CPD scaled conjugate cogradient.
    E = cache.factorizations.residual{f};
    offset = cache.factorizations.offset{f};
    grad = zeros(offset(end)-1,1);
    for n = 1:length(z)
        tmp = full(mtkrprod(E,z,n));
        grad(offset(n):offset(n+1)-1) = tmp(:);
    end
    
end

% Model: block term decomposition -----------------------------------------

function state_btd(~,~)
end

function fval = objfun_btd(f,z)
    
    % BTD objective function.
    T = model.factorizations{f}.data;
    isincomplete = cache.factorizations.isincomplete(f);
    issparse = cache.factorizations.issparse(f);
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
        T.val = fval;
        if ~isempty(T.matrix)
            T.matrix = sparse(double(T.sub{1}), ...
                double(1+idivide(T.ind-1,int64(size(T.matrix,1)))), ...
                double(fval),size(T.matrix,1),size(T.matrix,2));
        end
        cache.factorizations.residual{f} = T;
    else
        if issparse, size_tens = T.size;
        else size_tens = size(T); end
        cache.factorizations.residual{f} = reshape(fval,size_tens);
    end
    fval = 0.5*(fval(:)'*fval(:));
    
end

function grad = grad_btd(f,z)
    
    % BTD scaled conjugate cogradient.
    N = length(z{1})-1;
    E = cache.factorizations.residual{f};
    offset = cache.factorizations.offset{f};
    grad = zeros(offset(end)-1,1);
    cnt = 1;
    for r = 1:length(z)
        U = z{r}(1:N);
        S = conj(z{r}{end});
        for n = 1:N
            tmp = full(mtkronprod(E,U,n))* ...
                reshape(permute(S,[1:n-1 n+1:N n]),[],size(S,n));
            grad(offset(cnt):offset(cnt+1)-1) = tmp(:);
            cnt = cnt+1;
        end
        tmp = full(mtkronprod(E,U,0));
        grad(offset(cnt):offset(cnt+1)-1) = tmp;
        cnt = cnt+1;
    end
    
end

% Model: L2 regularization ------------------------------------------------

function state_regL2(f,firstrun)
    
    % Format right hand side, if available.
    if nargin == 2 && firstrun
         cache.factorizations.hasdata(f) = ...
             isfield(model.factorizations{f},'data');
         if cache.factorizations.hasdata(f) && ...
            ~iscell(model.factorizations{f}.data)
             model.factorizations{f}.data = ...
                 {model.factorizations{f}.data};
         end
         if cache.factorizations.hasdata(f)
             model.factorizations{f}.data = cellfun(@(d)full(d), ...
                 model.factorizations{f}.data,'UniformOutput',false);
         end
    end
    
end

function fval = objfun_regL2(f,z)
    % L2 regularization objective function.
    fval = 0;
    for i = 1:length(z)
        e = z{i}(:);
        if cache.factorizations.hasdata(f)
            e = e-model.factorizations{f}.data{i}(:);
        end
        fval = fval+(e'*e);
    end
    fval = 0.5*fval;
end

function grad = grad_regL2(f,z)
    % L2 regularization scaled conjugate cogradient.
    if cache.factorizations.hasdata(f)
        grad = cell2mat(cellfun(@(f,d)f(:)-d(:), ...
            z(:),model.factorizations{f}.data(:),'UniformOutput',false));
    else
        grad = cell2mat(cellfun(@(f)f(:),z(:),'UniformOutput',false));
    end
end

% Model: L1 regularization ------------------------------------------------

function state_regL1(f,firstrun)
    
    % Format right hand side, if available.
    if nargin == 2 && firstrun
         cache.factorizations.hasdata(f) = ...
             isfield(model.factorizations{f},'data');
         if cache.factorizations.hasdata(f) && ...
            ~iscell(model.factorizations{f}.data)
             model.factorizations{f}.data = ...
                 {model.factorizations{f}.data};
         end
         if cache.factorizations.hasdata(f)
             model.factorizations{f}.data = cellfun(@(d)full(d), ...
                 model.factorizations{f}.data,'UniformOutput',false);
         end
         if ~isfield(options,'mu')
             if cache.factorizations.hasdata(f)
                m = max(abs(serialize(model.factorizations{f}.data)));
             else
                m = 0;
             end
             % For x <= mu, use a smooth rational approximation of abs(x).
             % For x > mu, use abs(x).
             options.Mu = max(m,1)/100;
         end
    end
    
end

function fval = objfun_regL1(f,z)
    % Approximate L1 regularization objective function.
    fval = 0;
    for i = 1:length(z)
        e = z{i}(:);
        if cache.factorizations.hasdata(f)
            e = e-model.factorizations{f}.data{i}(:);
        end
        far = abs(e) > options.Mu;
        e(far) = abs(e(far));
        e2 = e(~far);
        e2 = e2.*conj(e2);
        e(~far) = 2*options.Mu*e2./(e2+options.Mu^2);
        fval = fval+sum(e);
    end
    fval = 0.5*fval;
end

function grad = grad_regL1(f,z)
    % Approximate L1 regularization scaled conjugate cogradient.
    if cache.factorizations.hasdata(f)
        grad = cell2mat(cellfun(@(f,d)f(:)-d(:), ...
            z(:),model.factorizations{f}.data(:),'UniformOutput',false));
    else
        grad = cell2mat(cellfun(@(f)f(:),z(:),'UniformOutput',false));
    end
    far = abs(grad) > options.Mu;
    grad(far) = grad(far)./(2*abs(grad(far)));
    grad(~far) = 2*options.Mu^3*grad(~far)./ ...
        (grad(~far).*conj(grad(~far))+options.Mu^2).^2;
end

end
