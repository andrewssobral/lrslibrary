function [U,S,sv] = mlsvd(T,size_core,perm)
%MLSVD (Truncated) multilinear singular value decomposition.
%   [U,S,sv] = mlsvd(T) computes the factor matrices U{1}, ..., U{N} and
%   all-unitary and ordered core tensor S belonging to an MLSVD of the N-th
%   order tensor T. The mode-n multilinear singular values are stored in
%   sv{n}. The vector sv{n} contains the Frobenius norms of the mode-n
%   slices of the core tensor S, e.g. S(:,:,k) are the mode-3 slices of S.
%
%   [U,S] = mlsvd(T,size_core) and mlsvd(T,size_core,perm) return a
%   sequentially truncated MLSVD, in which the core tensor S has dimensions
%   size_core and is no longer all-unitary in general. If the vector perm
%   is not specified or empty, the modes are processed in the order that
%   heuristically minimizes the computational complexity. Otherwise, they
%   are processed in the order defined by the vector perm, e.g. 1:ndims(T).
%
%   [U,S] = mlsvd(T,tol) and mlsvd(T,tol,perm) return a sequentially
%   truncated MLSVD, in which the core tensor has dimensions such that the
%   relative error in Frobenius norm, frob(T-lmlragen(U,S))/frob(T), does
%   not exceed the tolerance tol, chosen as a value in the interval [0,1].
%   The dimensions of the core tensor S are chosen heuristically, in some
%   cases better results may be obtained by choosing them manually.
%
%   [U,S] = mlsvd(T,size_core,0) and mlsvd(T,tol,0) return a (parallelly)
%   truncated MLSVD, defined by either the truncation dimensions size_core
%   or the relative error tolerance tol (see above).
%
%   See also lmlra, tmprod, svd.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be),
%            Nick Vannieuwenhoven (Nick.Vannieuwenhoven@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L.R. Tucker, "The Extension of Factor Analysis to Three-Dimensional
%       Matrices", in Contributions to mathematical psychology,
%       H. Gulliksen and  N. Frederiksen, eds., Holt, Rinehart & Winston,
%       NY, 1964, pp. 109-127.
%   [2] L.R. Tucker, "Some Mathematical Notes on Three-Mode Factor
%       Analysis", Psychometrika, Vol. 31, 1966, pp. 279-311.
%   [3] L. De Lathauwer, B. De Moor, J. Vandewalle, "A Multilinear Singular
%       Value Decomposition", SIAM J. Matrix Anal. Appl., Vol. 21, No. 4,
%       April 2000, pp. 1253-1278.
%   [4] N. Vannieuwenhoven, R. Vandebril, K. Meerbergen, "A new truncation
%       strategy for the higher-order singular value decomposition", SIAM
%       J. Sci. Comput., Vol. 34, No. 2, April 2012, pp. A1027-A1052.

% Check the tensor T.
N = ndims(T);
size_tens = size(T);

% Check the dimensions of the core tensor and permutation vector.
isSTMLSVD = false;
isTMLSVD = false;
if nargin >= 2
    isSTMLSVD = true;
    size_core = size_core(:)';
    if length(size_core) == 1
        tol = size_core;
        if tol < 0 || tol > 1
            error('mlsvd:tol','Tolerance should be in [0,1].');
        end
    elseif length(size_core) == N
        tol = nan;
        if any(size_core > size_tens)
            error('mlsvd:size_core','size_core should be <= size(T).');
        end
    else
        error('mlsvd:arg2','length(size_core) should be ndims(T).');
    end
end
if nargin == 3 && ~isempty(perm) && ~isstruct(perm)
    if length(perm) == 1 && ~perm
        isTMLSVD = true;
        isSTMLSVD = false;
    elseif length(perm) ~= N
        error('mlsvd:perm', ...
             ['The vector perm should contain all integers 1:N, or, ' ...
              'in the case of a T-MLSVD, be 0.']);
    end
end

% Choose a permutation order that heuristically minimizes the computational
% complexity (cf. [4]).
if ~isSTMLSVD
    % MLSVD or T-MLSVD.
    perm = 1:N;
elseif nargin < 3 || ~isnumeric(perm) || ...
       (isnumeric(perm) && length(perm) ~= N)
    % ST-MLSVD and perm not supplied.
    [~,perm] = sort(size_tens);
end

% Initialize the output.
U = cell(1,N);
sv = cell(1,N);
S = T;

% Process the modes in the order perm.
relerr = 0;
for n = 1:length(perm)
    
    % Compute SVD of the mode-perm(n) unfolding.
    [U{perm(n)},s,v] = svd(tens2mat(S,perm(n)),'econ');
    sv{perm(n)} = diag(s);
    
    % Truncate the left subspace if (S)T-MLSVD.
    if isSTMLSVD || isTMLSVD
        if isnan(tol)
            Un = U{perm(n)};
            U{perm(n)} = Un(:,1:min(size_core(perm(n)),size(Un,2)));
        else
            if n == 1, T2 = sum(sv{perm(n)}.^2); end
            cs = cumsum(flipud(sv{perm(n)}.^2));
            idx = find(relerr+sqrt(cs/T2) < tol/(N-n+1),1,'last');
            if ~isempty(idx)
                relerr = relerr+sqrt(cs(idx)/T2);
                Un = U{perm(n)};
                U{perm(n)} = Un(:,1:end-idx);
            end
        end
    end
    
    % If ST-MLSVD, use core tensor S instead of T for the next mode.
    size_core(perm(n)) = size(U{perm(n)},2);
    size_tens(perm(n)) = size_core(perm(n));
    if isSTMLSVD
        s = sv{perm(n)}.';
        v = bsxfun(@times, ...
                conj(v(:,1:size_core(perm(n)))), ...
                s(1:size_core(perm(n))));
        S = mat2tens(v,size_tens,[],perm(n));
    end
    
end

% Compute the core tensor if not already available.
if ~isSTMLSVD, S = tmprod(T,U,1:N,'H'); end
