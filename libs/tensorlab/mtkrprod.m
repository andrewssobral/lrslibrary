function M = mtkrprod(T,U,n,conjugate)
%MTKRPROD Compute a matricized tensor Khatri-Rao product.
%   mtkrprod(T,U,n) computes the product
%
%      tens2mat(T,n)*conj(kr(U([end:-1:n+1 n-1:-1:1])))
%
%   without permuting the tensor T. Note that for improved performance, it
%   is advisable for the largest two dimensions of T to be the first and
%   last modes of T.
%
%   mtkrprod(T,U,n,false) does not apply the complex conjugate to the
%   factor matrices U{n} before computing the matricized tensor Khatri-Rao
%   product.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] N. Vannieuwenhoven, N. Vanbaelen, K. Meerbergen, R. Vandebril,
%       "The dense multiple-vector tensor-vector product: An initial
%       study," Technical Report TW635, Dept. of Computer Science,
%       KU Leuven, 2013.
%   [2] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 4, conjugate = true; end
data = T;
size_tens = cellfun('size',U,1);
if isstruct(T)
    % Incomplete or sparse tensor.
    T = data.matrix;
    alloc = @sparse;
else
    % Full tensor.
    alloc = @zeros;
end
N = length(size_tens);
idx = [1:n-1 n+1:N];
R = size(U{idx(1)},2);

% Determine if large-scale version of the algorithm should be executed.
left = n == N || (n > 1 && size_tens(1) > size_tens(end));
if left, mem = 8*R*prod(size_tens(2:end));
else mem = 8*R*prod(size_tens(1:end-1));
end
largescale = mem > 2e9 || (isstruct(data) && isempty(data.matrix));
if largescale && ~isstruct(data)
    error('mtkrprod:mem',['Intermediate result is too large. Try ' ...
        'to format data with fmt first so that the large-scale ' ...
        'version of the algorithm can be applied.']);
end

if largescale
    
    % The large scale implementation fires if the required memory would
    % otherwise be larger than 2GB and assumes the dataset has been
    % formatted by fmt.
    idx = [1:n-1 n+1:N];
    if conjugate, U = cellfun(@conj,U,'UniformOutput',false); end
    if n == 0, M = zeros(1,R);
    else M = zeros(size_tens(n),R); end
    for r = 1:R
        s = data.val;
        for j = 1:length(idx)
            s = s.*U{idx(j)}(data.sub{idx(j)},r);
        end
        if n == 0, M(r) = sum(s);
        else M(:,r) = accumarray(double(data.sub{n}),s,[size_tens(n) 1]);
        end
    end
    
else

    if left

        % Apply tensor matrix product.
        if conjugate, M = U{1}'*reshape(T,size_tens(1),[]);
        else M = U{1}.'*reshape(T,size_tens(1),[]); end

        % Choose order of operations.
        perm = ones(1,N);
        l = 2; r = N;
        for i = 2:N
            if r == n || (l < n && size_tens(l) > size_tens(r))
                perm(i) = l; l = l+1;
            else
                perm(i) = r; r = r-1;
            end
        end

        % Apply structure-exploiting matricized tensor Khatri-Rao product.
        for i = 2:length(perm)
            mode = perm(i);
            if mode == n, continue; end
            if conjugate, Umode = conj(U{mode});
            else Umode = U{mode}; end
            tmp = M;
            M = alloc(R,size(M,2)/size_tens(mode));
            if mode < n
                for j = 1:R
                    tmp2 = Umode(:,j).'* ...
                        reshape(tmp(j,:),size_tens(mode),[]);
                    M(j,:) = tmp2(:);
                end
            elseif mode > n
                for j = 1:R
                    tmp2 = reshape(tmp(j,:),[],size_tens(mode))*Umode(:,j);
                    M(j,:) = tmp2(:);
                end
            end 
        end

        % Permute output.
        M = M.';

    else

        % Apply tensor matrix product.
        if conjugate, M = reshape(T,[],size_tens(end))*conj(U{end});
        else M = reshape(T,[],size_tens(end))*U{end}; end

        % Choose order of operations.
        perm = N*ones(1,N);
        l = 1; r = N-1;
        for i = 2:N
            if r == n || (l < n && size_tens(l) > size_tens(r))
                perm(i) = l; l = l+1;
            else
                perm(i) = r; r = r-1;
            end
        end

        % Apply structure-exploiting matricized tensor Khatri-Rao product.
        for i = 2:length(perm)
            mode = perm(i);
            if mode == n, continue; end
            if conjugate, Umode = conj(U{mode});
            else Umode = U{mode}; end
            tmp = M;
            M = alloc(size(M,1)/size_tens(mode),R);
            if mode < n
                for j = 1:R
                    tmp2 = Umode(:,j).'* ...
                        reshape(tmp(:,j),size_tens(mode),[]);
                    M(:,j) = tmp2(:);
                end
            elseif mode > n
                for j = 1:R
                    tmp2 = reshape(tmp(:,j),[],size_tens(mode))*Umode(:,j);
                    M(:,j) = tmp2(:);
                end
            end
        end

    end

end

end
