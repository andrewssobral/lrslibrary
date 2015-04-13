function [S,iperm] = tmprod(T,U,mode,transpose,saveperm)
%TMPROD Mode-n tensor-matrix product.
%   S = tmprod(T,U,mode) computes the tensor-matrix product of the tensor T
%   with the matrices U{1}, ..., U{N} along the modes mode(1), ...,
%   mode(N), respectively. Note that in this implementation, the vector
%   mode should contain distinct integers. The mode-n tensor-matrix
%   products are computed sequentially in a heuristically determined order.
%   A mode-n tensor-matrix product results in a new tensor S in which the
%   mode-n vectors of a given tensor T are premultiplied by a given matrix
%   U{n}, i.e., tens2mat(S,mode(n)) = U{n}*tens2mat(T,mode(n)).
%
%   S = tmprod(T,U,mode,'T') and tmprod(T,U,mode,'H') apply the mode-n
%   tensor-matrix product using the transposed matrices U{n}.' and
%   conjugate transposed matrices U{n}' respectively along mode(n).
%
%   [S,iperm] = tmprod(T,U,mode) and S = tmprod(T,U,mode,'saveperm') save
%   one permutation operation. In the former case, the tensor-matrix
%   product can then be recovered by permute(S,iperm).
%
%   See also tens2mat, mat2tens.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be),
%            Nick Vannieuwenhoven (Nick.Vannieuwenhoven@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check arguments.
if nargin < 4, transpose = 0; end
if nargin < 5, saveperm = false; end
if ischar(saveperm), saveperm = strcmpi(saveperm,'saveperm'); end
switch transpose
    case {'T','H'}, m = [2 1];
    otherwise, m = [1 2]; if nargin < 5, saveperm = transpose; end
end
if ~iscell(U), U = {U}; end
if length(U) ~= length(mode)
    error('tmprod:NumberOfProducts','length(U) should be length(mode).');
end
U = U(:)'; mode = mode(:)';
size_tens = ones(1,max(mode));
size_tens(1:ndims(T)) = size(T);
if any(cellfun('size',U,m(2)) ~= size_tens(mode))
    error('tmprod:U','size(T,mode(n)) should be size(U{n},%i).',m(2));
end

% Sort the order of the mode-n products.
[~,idx] = sort(size_tens(mode)./cellfun('size',U,m(1)));
mode = mode(idx);
U = U(idx);

% Compute the complement of the set of modes.
n = length(mode);
N = length(size_tens);
bits = ones(1,N);
bits(mode) = 0;
modec = 1:N;
modec = modec(logical(bits(modec)));

% Prepermute the tensor.
perm = [mode modec];
size_tens = size_tens(perm);
S = T; if any(mode ~= 1:n), S = permute(S,perm); end

% Cycle through the n-mode products.
for i = 1:n
    size_tens(1) = size(U{i},m(1));
    switch transpose
        case 'T'
            S = reshape(U{i}.'*reshape(S,size(S,1),[]),size_tens);
        case 'H'
            S = reshape(U{i}'*reshape(S,size(S,1),[]),size_tens);
        otherwise
            S = reshape(U{i}*reshape(S,size(S,1),[]),size_tens);
    end
    if i < n
        S = permute(S,[2:N 1]);
        size_tens = size_tens([2:N 1]);
    end
end

% Inverse permute the tensor, unless the user intends to do so himself.
iperm(perm([n:N 1:n-1])) = 1:N;
if nargout <= 1 && ~saveperm, S = permute(S,iperm); end
