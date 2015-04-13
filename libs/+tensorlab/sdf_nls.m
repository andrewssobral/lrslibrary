function [sol,output] = sdf_nls(model,options)
%SDF_NLS Structured data fusion by nonlinear least squares.
%   [sol,output] = sdf_nls(model) solves the data fusion problem described
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
%   sdf_nls(model,options) may be used to set the following options:
%
%      options.Algorithm =   - The desired optimization method.
%      [@nls_gncgs| ...
%       {@nls_gndl}|@nls_lm]
%      options.PC = true     - Whether or not to use a preconditioner, if
%                              available, for computing the Gauss-Newton
%                              step.
%      options.RelWeights =  - By supplying relative weights, the weights
%      ones(1,F)               options.Weights are computed as follows:
%                              options.Weights(f) = options.RelWeights(f)/
%                              (sum(options.RelWeights)*numel(data_f)) for
%                              each of the F factorizations in the model.
%      options.Weights       - The weight of the F factorizations in the
%                              data fusion model. The SDF objective
%                              function is \sum_f 0.5*options.Weights(f)*
%                              frob(model_f-data_f)^2, where model_f and
%                              data_f are the fth factorization and
%                              corresponding tensor. By default, weights
%                              are provided by options.RelWeights, but
%                              options.Weights has precedence if supplied.
%      options.<...>         - Parameters passed to the selected method,
%                              e.g., options.TolFun, options.TolX.
%                              See also help [options.Algorithm].
%
%   See also sdf_minf.

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
    funcs = {@nls_gndl,@nls_gncgs,@nls_lm};
    options.Algorithm = funcs{find(cellfun(xsfunc,funcs),1)};
end
if ~isfield(options,'MaxIter'), options.MaxIter = 5000; end
if ~isfield(options,'CGMaxIter'), options.CGMaxIter = 15; end
if ~isfield(options,'Display'), options.Display = 0; end
if ~isfield(options,'JHasFullRank'), options.JHasFullRank = false; end
if ~isfield(options,'PC'), options.PC = true; end
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
% these four functions.
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
            model.factorizations{I}.JHJx = @JHJx_cpd;
        case 'btd'
            model.factorizations{I}.state = @state_btd;
            model.factorizations{I}.objfun = @objfun_btd;
            model.factorizations{I}.grad = @grad_btd;
            model.factorizations{I}.JHJx = @JHJx_btd;
        case 'regL2'
            model.factorizations{I}.state = @state_regL2;
            model.factorizations{I}.objfun = @objfun_regL2;
            model.factorizations{I}.grad = @grad_regL2;
            model.factorizations{I}.JHJx = @JHJx_regL2;
        case 'regL1'
            model.factorizations{I}.state = @state_regL1;
            model.factorizations{I}.objfun = @objfun_regL1;
            model.factorizations{I}.grad = @grad_regL1;
            model.factorizations{I}.JHJx = @JHJx_regL1;
        otherwise
            error('sdf_nls:model','Model %s is not supported', ...
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
dF.JHF = @grad;
dF.JHJx = @JHJx;
if options.PC && isfield(options,'M') && isa(options.M,'function_handle')
    dF.M = options.M;
end
[z,output] = options.Algorithm(@objfun,dF,model.variables,options);
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

function x = derivexpand(r)
% Linearly expand variables using their transformations' Jacobians.

    % Save some references for speed.
    factors = model.factors;
    voffset = cache.variables.offset;
    isconst = cache.factors.isconst;
    sequence = cache.factors.sequence;
    state = cache.factors.state;
    structure = cache.factors.structure;

    % For each ith factor...
    x = cell(size(factors));
    for i = 1:length(factors)
        % For each jth subfactor...
        sub = cell(size(factors{i}));
        for j = 1:numel(factors{i})
            if isconst{i}(j)
                % Constant subfactor.
                sub{j} = zeros(size(sequence{i}{j}{1}));
            else
                % For each kth transformation...
                seq = factors{i}{j};
                ftr = sequence{i}{j};
                stt = state{i}{j};
                dim = structure{i}{j}{1};
                sub{j} = r(voffset(seq{1}):voffset(seq{1}+1)-1);
                if isnumeric(dim)
                    if ~isempty(dim), sub{j} = reshape(sub{j},dim); end
                else sub{j} = deserialize(sub{j},dim);
                end
                for k = 2:length(seq)
                    task = stt{k-1};
                    task.l = [];
                    task.r = sub{j};
                    sub{j} = seq{k}(ftr{k-1},task);
                end
            end

        end
        if numel(sub) == 1, x{i} = sub{1};
        else x{i} = cell2mat(sub);
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

function y = JHJx(~,x)
    
    % Expand the variables into factors.
    x = derivexpand(x);
    
    % Compute (J(z)'*J(z))*x after expansion.
    y = zeros(cache.variables.offset(end)-1,1);
    for i = find(options.Weights(:).' ~= 0)
        
        % Apply fast matrix-vector product.
        tmpa = getfactors(i,cache.factors.expanded);
        tmpb = getfactors(i,x);
        JHJx = model.factorizations{i}.JHJx(i,tmpa,tmpb);
        
        % Contract the factor matrices into variables.
        if i == 1
            y = options.Weights(i)*derivcontract(i,y,JHJx);
        else
            y = y+options.Weights(i)*derivcontract(i,y,JHJx);
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

function state_cpd(f,firstrun)
% Can read from model, cache and options, can save state by writing to the
% structure cache.factorizations.state{f}.

    if nargin == 2 && firstrun
        
        % Store the fraction of known elements.
        if cache.factorizations.isincomplete(f)
            cache.factorizations.scale{f} = ...
                length(model.factorizations{f}.data.val)./...
                prod(model.factorizations{f}.data.size);
        end
        
        % Set Block-Jacobi preconditioner if
        % - the model is a single CPD and
        % - each factor consists of exactly one subfactor and
        % - no subfactor is transformed and
        % - every subfactor is a reference to a variable and
        % - all variable references are unique and
        % - all factor references are unique and
        % - an NLS algorithm is used.
        options.BJ = false;
        var = cellfun(@(v)v{1}{1},model.factors,'UniformOutput',false);
        ftr = cell2mat(model.factorizations{f}.factors);
        if options.PC && numel(model.factorizations) == 1 && ...
                all(cellfun(@(v)length(v) == 1,model.factors)) && ...
                all(cellfun(@(v)length(v{1}) == 1,model.factors)) && ...
                all(cellfun(@isscalar,var)) && ...
                length(unique(cell2mat(var))) == length(cell2mat(var))&&...
                length(unique(ftr)) == length(ftr) && ...
                ~isempty(strfind(func2str(options.Algorithm),'nls'))
            options.M = @pc_cpd;
            options.BJ = true;
        end
        
    end

    % Cache the factor matrices' Gramians.
    idx = cache.factorizations.serialized{f};
    N = length(idx);
    R = size(cache.factors.expanded{idx(1)},2);
    cache.factorizations.state{f}.UHU = zeros(N,R*R);
    for n = 1:N
        tmp = cache.factors.expanded{idx(n)};
        tmp = conj(tmp'*tmp);
        cache.factorizations.state{f}.UHU(n,:) = tmp(:);
    end
    
    % Optionally cache the inverses of the Gramians for the preconditioner.
    % In a faster language, this should be the Cholesky factor instead.
    if options.BJ
        cache.factorizations.state{f}.invW = cell(1,N);
        for n = 1:N
            tmp = cache.factorizations.state{f}.UHU([1:n-1 n+1:N],:);
            if N > 2, tmp = prod(tmp,1); end
            cache.factorizations.state{f}.invW{n} = ...
                inv(reshape(tmp,[R R]));
        end
    end
    
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

function y = JHJx_cpd(f,z,x)
    
    % CPD fast Jacobian's Gramian vector product.
    % Ignores the fact that the tensor might be incomplete.
    R = size(z{1},2);
    N = length(z);
    offset = cache.factorizations.offset{f};
    UHU = cache.factorizations.state{f}.UHU;
    XHU = zeros(R,R,N);
    y = zeros(offset(end)-1,1);
    for n = 1:N
        Wn = UHU([1:n-1 n+1:N],:);
        if N > 2, Wn = prod(Wn,1); end
        XHU(:,:,n) = conj(x{n}'*z{n});
        y(offset(n):offset(n+1)-1) = x{n}*reshape(Wn,[R R]);
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
    if cache.factorizations.isincomplete(f)
        y = y*cache.factorizations.scale{f};
    end
    
end

function x = pc_cpd(~,b)

    % Solve M*x = b, where M is a block-diagonal approximation for JHJ.
    % Equivalent to simultaneous ALS updates for each of the factors.
    x = zeros(size(b));
    offset = cache.factorizations.offset{1};
    invW = cache.factorizations.state{1}.invW;
    for n = 1:length(offset)-1
        idx = offset(n):offset(n+1)-1;
        tmp = reshape(b(idx),[],size(invW{1},1))*invW{n};
        x(idx) = tmp(:);
    end
    x = x/options.Weights;
    
    % If incomplete, approximate the effect of missing entries.
    if cache.factorizations.isincomplete(1)
        x = x/cache.factorizations.scale{1};
    end
    
end

% Model: block term decomposition -----------------------------------------

function state_btd(f,firstrun)
% Can read from model, cache and options, can save state by writing to the
% structure cache.factorizations.state{f}.
    
    if nargin == 2 && firstrun
        
        % Store the fraction of known elements.
        if cache.factorizations.isincomplete(f)
            cache.factorizations.scale{f} = ...
                length(model.factorizations{f}.data.val)./...
                prod(model.factorizations{f}.data.size);
        end
        
        % Set Block-Jacobi preconditioner if
        % - the model is a single BTD and
        % - each factor consists of exactly one subfactor and
        % - no subfactor is transformed and
        % - every subfactor is a reference to a variable and
        % - all variable references are unique and
        % - all factor references are unique and
        % - an NLS algorithm is used.
        options.BJ = false;
        var = cellfun(@(v)v{1}{1},model.factors,'UniformOutput',false);
        ftr = cache.factorizations.serialized{1};
        if options.PC && numel(model.factorizations) == 1 && ...
                all(cellfun(@(v)length(v) == 1,model.factors)) && ...
                all(cellfun(@(v)length(v{1}) == 1,model.factors)) && ...
                all(cellfun(@isscalar,var)) && ...
                length(unique(cell2mat(var))) == length(cell2mat(var))&&...
                length(unique(ftr)) == length(ftr) && ...
                ~isempty(strfind(func2str(options.Algorithm),'nls'))
            options.M = @pc_btd;
            options.BJ = true;
        end
        
    end
    
    % Cache the factor matrices' Gramians.
    U = getfactors(f,cache.factors.expanded);
    R = length(U);
    N = length(U{1})-1;
    [idx,jdx,kdx] = ndgrid(1:R,1:N,1:R);
    cache.factorizations.state{f}.UHU = ...
        arrayfun(@(i,n,j)U{i}{n}'*U{j}{n}, ...
        idx,jdx,kdx,'UniformOutput',false);
    
    % Optionally cache some results for the block-Jacobi preconditioner.
    if options.BJ
        [idx,jdx] = ndgrid(1:R,1:N);
        UHU = cache.factorizations.state{f}.UHU;
        cache.factorizations.state{f}.invSKS = arrayfun( ...
            @(r,n)inv(mtkronprod(U{r}{end},UHU(r,:,r),n)* ...
            conj(reshape(permute(U{r}{end},[1:n-1 n+1:N n]), ...
            [],size(U{r}{end},n)))),idx,jdx,'UniformOutput',false);
        cache.factorizations.state{f}.invUHU = arrayfun( ...
            @(r,n)inv(UHU{r,n,r}),idx,jdx,'UniformOutput',false);
    end

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

function y = JHJx_btd(f,z,x)
    
    % BTD fast Jacobian's Gramian vector product.
    % Ignores the fact that the tensor might be incomplete.
    R = length(z);
    N = length(z{1})-1;
    offset = cache.factorizations.offset{f};
    UHU = cache.factorizations.state{f}.UHU;
    [idx,jdx,kdx] = ndgrid(1:R,1:N,1:R);
    XHU = arrayfun(@(i,n,j)x{i}{n}'*z{j}{n}, ...
        idx,jdx,kdx,'UniformOutput',false);
    y = zeros(offset(end)-1,1);
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
    if cache.factorizations.isincomplete(f)
        y = y*cache.factorizations.scale{f};
    end
    
end

function x = pc_btd(~,b)

    % Solve M*x = b, where M is a block-diagonal approximation for JHJ.
    x = zeros(size(b));
    R = length(model.factorizations{1}.factors);
    N = length(model.factorizations{1}.factors{1})-1;
    offset = cache.factorizations.offset{1};
    invSKS = cache.factorizations.state{1}.invSKS;
    invUHU = cache.factorizations.state{1}.invUHU;
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
    x = x/options.Weights;
    
    % If incomplete, approximate the effect of missing entries.
    if cache.factorizations.isincomplete(1)
        x = x/cache.factorizations.scale{1};
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

function y = JHJx_regL2(~,~,x)
    % L2 regularization Hessian vector product.
    y = cell2mat(cellfun(@(f)f(:),x(:),'UniformOutput',false));
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

function y = JHJx_regL1(f,z,x)
    % Approximate L1 regularization Hessian vector product.
    if cache.factorizations.hasdata(f)
        hess = cell2mat(cellfun(@(f,d)f(:)-d(:), ...
            z(:),model.factorizations{f}.data(:),'UniformOutput',false));
    else
        hess = cell2mat(cellfun(@(f)f(:),z(:),'UniformOutput',false));
    end
    if ~isreal(hess)
        error('sdf_nls:regL1',['Please use sdf_minf when applying L1 ' ...
            'regularization on complex factors.']);
    end
    far = abs(hess) > options.Mu;
    hess(far) = 0;
    if any(hess)
        hess(~far) = 2*options.Mu^3*(options.Mu^2-3*hess(~far).^2)./ ...
            (options.Mu^2+hess(~far).^2).^3;
        y = hess.*cell2mat(cellfun(@(f)f(:),x(:),'UniformOutput',false));
    else
        y = hess;
    end
end

end
