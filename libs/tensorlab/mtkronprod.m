function M = mtkronprod(T,U,n,transpose)
%MTKRONPROD Compute a matricized tensor Kronecker product.
%   mtkronprod(T,U,n) computes the product
%
%      tens2mat(T,n)*conj(kron(U([end:-1:n+1 n-1:-1:1])))
%
%   without permuting the tensor T. Note that for improved performance, it
%   is advisable for the largest two dimensions of T to be the first and
%   last modes of T.
%
%   mtkronprod(T,U,n,'T') and mtkronprod(T,U,n,'H') transpose or complex
%   conjugate transpose the matrices U{n} before computing the matricized
%   tensor Kronecker product.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 4, transpose = 0; end
if ischar(transpose)
    if strcmp(transpose,'T'), transpose = 1;
    elseif strcmp(transpose,'H'), transpose = 1i; end
end
U = U(:).';
if transpose == 0, size_tens = cellfun('size',U,1);
else size_tens = cellfun('size',U,2); end
data = T;
if isstruct(T)
    % Incomplete or sparse tensor.
    T = data.matrix;
    alloc = @sparse;
else
    % Full tensor.
    alloc = @zeros;
end
N = length(size_tens);
if transpose == 0, size_core = cellfun('size',U,2);
else size_core = cellfun('size',U,1); end
if n > 0, size_core(n) = size_tens(n); end
R = prod(size_core([1:n-1 n+1:N]));

% Determine if large-scale version of the algorithm should be executed.
ratio = size_core./size_tens;
perm = zeros(1,N);
l = 1; r = N;
for i = 1:N
    if r == n || (l < n && ratio(l) < ratio(r)), perm(i) = l; l = l+1;
    else perm(i) = r; r = r-1;
    end
end
if n == 0, mem = 8*prod(size_tens);
else mem = 8*prod(size_tens([1:perm(1)-1 perm(1)+1:N]))*size_core(perm(1));
end
largescale = mem > 2e9 || (isstruct(data) && isempty(data.matrix));
if largescale && ~isstruct(data)
    error('mtkronprod:mem',['Intermediate result is too large. Try ' ...
        'to format data with fmt first so that the large-scale ' ...
        'version of the algorithm can be applied.']);
end

if largescale
    
    % The large scale implementation fires if the required memory would
    % otherwise be larger than 2GB and assumes the dataset has been
    % formatted by fmt.
    idx = [1:n-1 n+1:N];
    jdx = cell(1,N);
    switch transpose
        case 0,  U = cellfun(@(u)conj(u),U,'UniformOutput',false);
        case 1,  U = cellfun(@(u)u',U,'UniformOutput',false);
        case 1i, U = cellfun(@(u)u.',U,'UniformOutput',false);
    end
    if n == 0, M = zeros(1,R);
    else M = zeros(size_tens(n),R); end
    for r = 1:R
        s = data.val;
        [jdx{:}] = ind2sub(size_core([1:n-1 n+1:N]),r);
        for j = 1:length(idx)
            s = s.*U{idx(j)}(data.sub{idx(j)},jdx{j});
        end
        if n == 0, M(r) = sum(s);
        else M(:,r) = accumarray(double(data.sub{n}),s,[size_tens(n) 1]);
        end
    end
    
else

    % Apply structure-exploiting matricized tensor Kronecker product.
    M = T;
    cpl = cumprod([1 size_core(perm(perm < n))]); l = 1;
    cpr = cumprod([1 size_core(perm(perm > n))]); r = 1;
    for i = 1:length(perm)

        mode = perm(i);
        if mode < n

            tmp = reshape(M,cpl(l)*size_tens(mode),[]);
            M = alloc(cpl(l),size(tmp,2)*size_core(mode));
            for j = 1:cpl(l)
                idx = j:cpl(l):size(tmp,1);
                switch transpose
                    case 0,  tmp2 = U{mode}'*tmp(idx,:);
                    case 1,  tmp2 = conj(U{mode})*tmp(idx,:);
                    case 1i, tmp2 = U{mode}*tmp(idx,:);
                end
                M(j,:) = tmp2(:);
            end
            l = l+1;

        elseif mode > n

            tmp = reshape(M,[],cpr(r)*size_tens(mode));
            M = alloc(size(tmp,1)*size_core(mode),cpr(r));
            for j = 1:cpr(r)
                idx = (1:size_tens(mode))+(j-1)*size_tens(mode);
                switch transpose
                    case 0,  tmp2 = tmp(:,idx)*conj(U{mode});
                    case 1,  tmp2 = tmp(:,idx)*U{mode}';
                    case 1i, tmp2 = tmp(:,idx)*U{mode}.';
                end
                M(:,j) = tmp2(:);
            end
            r = r+1;

        end

    end

    % Permute and reshape output.
    if n > 1
        if issparse(M)
            idx = find(M);
            [i,j,k] = ind2sub([cpl(end) size_tens(n) cpr(end)],idx);
            ik = sub2ind([cpl(end) cpr(end)],i,k);
            M = sparse(j,ik,full(M(idx)),size_tens(n),cpl(end)*cpr(end));
        else
            M = reshape(M,[cpl(end) size_tens(n) cpr(end)]);
            M = reshape(permute(M,[2 1 3]),size_tens(n),[]);
        end
    else
        if n == 0, M = reshape(M,1,[]);
        else M = reshape(M,size_tens(n),[]); end
    end

end

end
