function C = scov(X,shift,prewhiten)
%SCOV Shifted covariance matrices.
%   C = scov(X,shift) computes shifted covariance matrices C(:,:,k) of a 
%   matrix X in which each row is an observation at a certain timestep and
%   each column is a variable. Herein,
%
%      C(i,j,k) = E[xi(t).*conj(xj(t+shift(k))]
%
%   where the expectation E is approximated by a mean that is normalized by
%   size(X,1)-1-shift(k), xi(t) is the i-th mean centered variable,
%   X(:,i)-mean(X(:,i)), and xj(t+shift(k)) is the j-th mean centered
%   variable shifted by shift(k) timesteps. The vector shift should contain
%   integers between 0 and size(X,1)-2.
%
%   C = scov(X) or C = scov(X,[]) chooses as many shifts as there are
%   variables, if there are enough observations, otherwise length(shifts)
%   is equal to size(X,1)-2. The shifts are spaced one timelag apart.
%
%   C = scov(X,shift,'prewhiten') applies a linear transformation to the
%   columns of X so that the covariance matrix of the new matrix (without
%   shift) is the identity matrix before computing its shifted covariance
%   matrices.
%
%   See also cov, cum4.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] A. Belouchrani, K. Abed-Meraim, J.-F. Cardoso, E. Moulines, "A
%       Blind Source Separation Technique Using Second Order Statistics",
%       IEEE Trans. Sig. Proc., Vol. 45, No. 2, February 1997, pp. 434-–444.
%   [2] Lieven De Lathauwer, Joséphine Castaing, "Blind Identification of
%       Underdetermined Mixtures by Simultaneous Matrix Diagonalization",
%       IEEE Trans. Sig. Proc., Vol. 56, No. 3, March 2008, pp. 1096-1105.

% Store the dimensions of the input matrix.
[n,v] = size(X);

% Check the shifts and prewhiten option.
if nargin < 2 || isempty(shift)
    shift = 0:min(v-1,n-2);
else
    if ~any(size(shift) == length(shift)) || ...
        any(shift < 0) || any(shift > n-2)
        error('scov:shift',['The argument shift should be a vector of ' ...
                            'integers between 0 and size(X,1)-2.']);
    end
end
if nargin < 3, prewhiten = false; end
if ischar(prewhiten), prewhiten = strcmpi(prewhiten,'prewhiten'); end

% Center the variables.
X = bsxfun(@minus,X,mean(X));

% Apply a prewhitening to X if requested.
if prewhiten
    [U,S,~] = svd(X,'econ');
    X = U*(S*pinv(S))*sqrt(n-1);
end

% Compute shifted variance-covariance matrices.
C = zeros(size(X,2),size(X,2),length(shift));
for i = 1:length(shift)
    t = shift(i);
	C(:,:,i) = conj(X(1:(end-t),:)'*X((t+1):end,:))/(n-1-t);
end
