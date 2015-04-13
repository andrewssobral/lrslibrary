function r = mlrank(T,tol)
%MLRANK Multilinear rank.
%   mlrank(A,tol) provides an estimate of the number of linearly
%   independent mode-n vectors in the tensor T, for n = 1, ..., ndims(T).
%   More specifically, r(n) is the number of mode-n multilinear singular
%   values that are greater than tol(n)*max(sv{n}), where sv{n} are the
%   mode-n multilinear singular values of T. If tol is a scalar, it is set
%   equal to tol*ones(1,ndims(T)).
%
%   mlrank(A) uses the tolerance tol(n) = max(size(tens2mat(T,n)))*eps.
%
%   See also mlsvd, rank.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Incomplete/sparse tensors not supported yet.
T = fmt(T,true);
if isstruct(T)
    error('This method does not yet support incomplete/sparse tensors.');
end

% Compute the multilinear singular values of T.
[~,~,sv] = mlsvd(T);

% Set the tolerance.
if nargin < 2, tol = max(size(T),numel(T)./size(T))*eps(class(T)); end
tol = tol(:).';
if length(tol) == 1, tol = tol*ones(1,ndims(T)); end
if length(tol) ~= ndims(T)
    error('mlrank:tol', ...
          'tol should be a scalar, or a vector of length ndims(T).');
end

% Compute the multilinear rank.
r = arrayfun(@(n)sum(sv{n} > tol(n)*max(sv{n})),1:ndims(T));
