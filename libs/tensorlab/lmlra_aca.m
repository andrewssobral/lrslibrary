function [U,S,output] = lmlra_aca(T,size_core,options)
%LMLRA_ACA LMLRA by adaptive cross-approximation.
%   [U,S] = lmlra_aca(T,size_core) computes factor matrices U{1},
%   ..., U{N} of dimensions size(T,n)-by-size_core(n) and a core tensor of
%   dimensions size_core by adaptive cross approximation.
%
%   lmlra_aca(T,size_core,options) may be used to set the following
%   options:
%
%      options.Display = false - Set to true to enable printing output
%                                information to the command line.
%      options.TolSV = 1e-12   - The relative singular value tolerance in
%                                determining the factor matrices.
%
%   See also lmlragen.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] Cesar F. Caiafa, Andrzej Cichocki, "Generalizing the column-row
%       matrix decomposition to multi-way arrays," Linear Algebra and its
%       Applications, Vol. 433, No. 3, 2010, pp. 557-573.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

% Check the options structure.
if nargin < 3, options = struct; end
if ~isfield(options,'Display'), options.Display = false; end
if ~isfield(options,'MaxTries'), options.MaxTries = 3; end
if ~isfield(options,'TolSV'), options.TolSV = 1e-12; end

% Preprocess data.
T = fmt(T,true);
N = length(size_core);
idx = cell(1,N);
if isstruct(T), size_tens = T.size;
else size_tens = [size(T) ones(1,length(size_core)-ndims(T))];
end

% Get primary indices.
if isstruct(T)
    
    % Count how many known elements per mode-n slice.
    cnt = arrayfun(@(n)accumarray(T.sub{n},1,[size_tens(n) 1]),1:N, ...
        'UniformOutput',false);
    
    % Take indices (i1,...,iN) corresponding to largest geometric mean
    % between fraction of known elements and the ratio of the absolute
    % value of the element w.r.t. the largest element in absolute value.
    score = cnt{1}(T.sub{1});
    for n = 2:N, score = score+cnt{n}(T.sub{n}); end
    score = score./sum(prod(size_tens)./size_tens);
    [~,tmp] = max(abs(T.val).*score);
    idx = arrayfun(@(n)double(T.sub{n}(tmp)),1:N,'UniformOutput',false);
    
    % Get fill-in values.
    fill = arrayfun(@(n)accumarray(T.sub{n},T.val,[size_tens(n) 1])./ ...
        cnt{n},1:N,'UniformOutput',false);
    
else
    
    % Take indices (i1,...,iN) corresponding to the largest element in T.
    for n = 1:N
        score = reshape(abs(T),[prod([1 size_tens(1:n-1)]) ...
            size_tens(n) prod([1 size_tens(n+1:N)])]);
        [~,idx{n}] = max(max(max(score,[],1),[],3));
    end
    
end

% Initialize adaptive cross approximation.
% Keep track of core tensor S and large factor matrices C, which together
% define the factor matrices U.
tries = 1;
cross = idx;
C = arrayfun(@(n)zeros(zeros(1,N)),1:N,'UniformOutput',false);
S = zeros(zeros(1,N));
update();
if options.Display
    if isstruct(T), val = abs(T.val(tmp)); else val = abs(T(idx{:})); end
    fprintf('size(S) = [%i',length(idx{1}));
    fprintf(' %i',cellfun(@length,idx(2:end)));
    fprintf('] (abs(T(%i',cross{1});
    fprintf(',%i',cell2mat(cross(2:end)));
    fprintf(')) = %g)\n',val);
end

% Modified version of the adaptive cross approximation algorithm of [1].
while ~all(cellfun(@length,idx) == size_core)
    
    % Compute relative errors in currently selected cross.
    err = nan(1,N);
    obj = nan(1,N);
    loc = nan(1,N);
    for n = 1:N
        if length(idx{n}) == size_core(n) || size_tens(n) == 1
            continue;
        end
        res = abs(residual(cross,n));
        res(idx{n}) = 0;
        if isstruct(T)
            [obj(n),loc(n)] = max(res(:).*cnt{n}*size_tens(n));
        else
            [obj(n),loc(n)] = max(res);
        end
        err(n) = res(loc(n));
    end
    
    % Place a new cross at the location of the largest relative error.
    [~,opt] = max(obj);
    if err(opt) > 100*eps
        idx{opt}(end+1) = loc(opt);
        update();
        tries = 1;
        cross = cellfun(@(i)i(end),idx,'UniformOutput',false);
        if options.Display
            fprintf('size(S) = [%i',length(idx{1}));
            fprintf(' %i',cellfun(@length,idx(2:end)));
            if numel(S) == 2
                fprintf('] (abs(T(%i',cross{1});
                fprintf(',%i',cell2mat(cross(2:end)));
                fprintf(')) = %g)\n',err(opt));
            else
                fprintf('] (Cross @ [%i',cross{1});
                fprintf(' %i',cell2mat(cross(2:end)));
                fprintf('] with rel.err. = %g)\n',err(opt));
            end
        end
    elseif isstruct(T) && tries <= options.MaxTries
        tries = tries+1;
        cross = cellfun(@(i)i(randi(length(i),1)),idx,'UniformOutput',0);
        if options.Display
            fprintf('size(S) = [%i',length(idx{1}));
            fprintf(' %i',cellfun(@length,idx(2:end)));
            fprintf('] (Trying random cross @ [%i',cross{1});
            fprintf(',%i',cell2mat(cross(2:end)));
            fprintf('])\n');
        end
    else
        break;
    end
    
end

% Build LMLRA.
U = build();
output.idx = idx;

function n = update()
    n = find(arrayfun(@(n)size(S,n)~=length(idx{n}),1:N),1,'last');
    jdx = idx;
    jdx{n} = idx{n}(end);
    pad = subsref(jdx,true,false);
    pad(~isfinite(pad)) = 1;
    kdx = repmat({':'},1,N);
    kdx{n} = size(S,n)+1;
    S(kdx{:}) = pad;
    for m = 1:N
        if m == n && numel(S) > 1, continue; end
        jdx = idx;
        jdx{n} = idx{n}(end);
        jdx{m} = 1:size_tens(m);
        pad = subsref(jdx,true,true);
        pad(~isfinite(pad)) = 1;
        kdx = repmat({':'},1,N);
        if numel(S) > 1, kdx{n} = size(C{m},n)+1; end
        C{m}(kdx{:}) = pad;
    end
end

function U = build(n)
    U = cell(1,N);
    for m = 1:N
        jdx = repmat({':'},1,N);
        if nargin == 1 && m ~= n, jdx{m} = idx{m}(end); end
        [u,s,v] = svd(tens2mat(S,m),'econ');
        s = diag(s);
        r = sum(s/s(1) > options.TolSV);
        U{m} = (tens2mat(C{m}(jdx{:}),m)*v(:,1:r))* ...
            bsxfun(@rdivide,u(:,1:r)',s(1:r));
        if nargin == 1 && m ~= n, U{m} = U{m}(:); end
    end
end

function res = residual(idx,n)
    idx{n} = 1:size_tens(n);
    res = subsref(idx);
    if numel(S) > 1
        u = build(n);
        v = res(:);
        res = v-u{n}*mtkrprod(S,u,n,false);
        idx = abs(v) > 1;
        res(idx) = res(idx)./abs(v(idx));
    end
end

function sub = subsref(idx,fil,sgn)
    if nargin == 1, fil = false; end
    if isstruct(T)
        sub = nan(cellfun(@length,idx));
        idx = cellfun(@int64,idx,'UniformOutput',false);
        [idx{:}] = ndgrid(idx{:});
        ind = sub2ind(size_tens,idx{:});
        indismem = ismember(T.ind,ind);
        Tind = T.ind(indismem);
        Tval = T.val(indismem);
        subismem = ismember(ind,Tind);
        [~,idx] = sort(ind(subismem));
        jdx(idx) = 1:sum(subismem(:));
        sub(subismem) = Tval(jdx);
        if fil
            idx = cell(1,length(size_tens));
            [idx{:}] = ind2sub(size_tens,find(~subismem));
            filled = fill{1}(idx{1});
            for m = 2:N, filled = filled.*fill{m}(idx{m}); end
            if isreal(T.val)
                if sgn, sgn = sign(filled); end
                filled = nthroot(abs(filled),N);
                if isnumeric(sgn), filled = sgn.*filled; end
            else
                filled = filled.^(1/N);
            end
            sub(~subismem) = filled;
        end
    else
        sub = T(idx{:});
    end
end

end
