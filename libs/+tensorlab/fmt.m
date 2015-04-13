function T = fmt(T,internal)
%FMT Format data set.
%   T = fmt(T) converts the tensor T into the format used by Tensorlab,
%   exploiting sparsity or incompleteness where possible. There are two
%   cases:
%
%   1. T is an array. If T contains many zeros or if it contains one or
%   more NaNs, it is converted to a structure (see below).
%
%   2. T is a structure. The following fields, if missing, will be
%   completed or formatted properly:
%
%      T.size       - The size of the tensor.
%      T.val        - The vector of nonzero or known elements.
%      T.ind        - The linear indices of the nonzero or known elements.
%      T.sub        - Either an array of subindices of the nonzero or known
%                     elements with one row per nonzero or known element.
%                     Or, a cell array with one vector of subindices per
%                     mode in T.
%      T.matrix     - A sparse matrix representation of the mode-1
%                     unfolding of T.
%      T.incomplete - Boolean which sets if the tensor T is incomplete.
%      T.sparse     - Boolean which sets if the tensor T is sparse.
%
%   See also sdf_minf, sdf_nls.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

% Called internally?
if nargin < 2, internal = false; end

% Check if dense matrix should be converted to incomplete/sparse format.
tic;
if isnumeric(T)

    size_tens = size(T);
    finite = isfinite(T);
    if ~all(finite(:))
        
        % Incomplete format.
        matrix = T;
        matrix(~finite) = 0;
        matrix = sparse(reshape(matrix,size_tens(1),[]));
        T = struct('matrix',matrix,'size',size_tens, ...
            'ind',find(finite),'val',T(finite), ...
            'incomplete',true,'sparse',false);

    elseif nnz(T)/numel(T) < 0.05
        
        % Sparse format.
        matrix = sparse(reshape(T,size_tens(1),[]));
        T = struct('matrix',matrix,'size',size_tens, ...
            'ind',find(T),'val',reshape(full(T(T~=0)),[],1),...
            'incomplete',false,'sparse',true);

    end
    
end

% If the T is a structure, fill in the blanks.
if isstruct(T)
    
    % Incomplete/sparse T must have a size field.
    if ~isfield(T,'size')
        error('fmt:size','Data is missing field ''size''.');
    end
    
    % If there is a sub field, format it in ind2sub style.
    if isfield(T,'sub') && isnumeric(T.sub)
        T.sub = mat2cell(T.sub, ...
            size(T.sub,1),ones(1,size(T.sub,2)));
    end
    
    % Fill in the ind, sub and val field.
    if isfield(T,'ind') && ~isfield(T,'sub')
        if ~isa(T.ind,'double'), T.ind = double(T.ind); end
        T.sub = cell(1,length(T.size));
        [T.sub{:}] = ind2sub(T.size,T.ind);
    end
    if isfield(T,'sub') && ~isfield(T,'ind')
        if ~all(cellfun(@(s)isa(s,'double'),T.sub))
            T.sub = cellfun(@(s)double(s),T.sub,'UniformOutput',false);
        end
        T.ind = sub2ind(T.size,T.sub{:});
    end
    if isfield(T,'matrix')
        if ~isfield(T,'ind')
            T.ind = find(T.matrix);
        end
        if ~isfield(T,'sub')
            T.sub = cell(1,length(T.size));
            [T.sub{:}] = ind2sub(T.size,T.ind);
        end
        if ~isfield(T,'val')
            T.val = full(T.matrix(T.ind));
        end
    end
    
    % Store values as double for compatibility.
    if ~isa(T.val,'double') && ~isa(T.val,'single')
        T.val = double(T.val);
    end
    
    % Fill in the matrix field as a mode-1 unfolding of the T.
    % This way, sparse matrix operations can be exploited.
    if ~isfield(T,'matrix')
        if ~isfield(T,'val')
            error('fmt:val','Data is missing field ''val''.');
        end
        [~,maxsize] = computer();
        if 1e6*prod(T.size)/min(T.size) > maxsize
            warning('fmt:maxsize',['Tensor is too large to ' ...
                'be stored using sparse matrices. This may lead ' ...
                'to compatibility issues.']);
            T.matrix = [];
        else
            if ~all(cellfun(@(s)isa(s,'double'),T.sub))
                T.sub = cellfun(@(s)double(s),T.sub,'UniformOutput',false);
            end
            row = T.sub{1};
            col = sub2ind([T.size(2:end) 1],T.sub{2:end});
            T.matrix = sparse(row,col,double(T.val), ...
                T.size(1),prod(T.size(2:end)));
        end
    end
    
    % Store linear indices as int64 for convenience.
    if ~isa(T.ind,'int64')
        T.ind = int64(T.ind);
    end
    
    % Store subindices as int32 to save space.
    if ~all(cellfun(@(s)isa(s,'int32'),T.sub))
        T.sub = cellfun(@(s)int32(s),T.sub,'UniformOutput',false);
    end
    
    % Sort entries for increased performance.
    if ~internal && ~issorted(T.ind)
        [T.ind,idx] = sort(T.ind);
        T.val = T.val(idx);
        T.sub = cellfun(@(s)s(idx),T.sub,'UniformOutput',false);
    end
    
    % Determine type of T.
    if isfield(T,'incomplete') && isfield(T,'sparse') && ...
        T.incomplete && T.sparse
        error('fmt:type','Tensor can not be both incomplete and sparse.');
    elseif isfield(T,'sparse') && ~isfield(T,'incomplete')
        T.incomplete = ~T.sparse;
    elseif isfield(T,'incomplete') && ~isfield(T,'sparse')
        T.sparse = ~T.incomplete;
    elseif ~isfield(T,'sparse') && ~isfield(T,'incomplete')
        warning('fmt:type', ...
            'Assuming T is incomplete (instead of sparse).');
        T.incomplete = true;
        T.sparse = false;
    end
    
end

% Show warning if conversion time took too long.
t = toc;
if internal && t > 0.3
    warning('sdf:fmt',['It is recommended to format this tensor ' ...
                'with fmt (conversion took %gs).'],t);
end
