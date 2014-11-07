function [Y,N] = noisy(X,SNR)
%NOISY Generate a noisy version of a given array.
%   [Y,N] = noisy(X,SNR) computes a noisy version of X as Y = X + N, where
%   the noise term N is generated as sigma*randn(size(X)) if X is real and
%   as sigma*(randn(size(X))+randn(size(X))*1i) if X is complex. The scalar
%   sigma is chosen such that 10*log10((X(:)'*X(:))/(N(:)'*N(:))) = SNR dB.
%   By default, SNR is 20. If X is a cell array, a noisy version of each of
%   its elements is recursively computed and returned in the cell array Y.
%
%   See also randn, rand.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check the options structure.
if nargin < 2, SNR = 20; end

% Add noise recursively.
[Y,N] = addnoise(X);

function [X,N] = addnoise(X)
    if iscell(X)
        N = X;
        for i = 1:length(X)
            [X{i},N{i}] = addnoise(X{i});
        end
        return;
    end
    if isstruct(X), size_x = size(X.val); real = isreal(X.val);
    else size_x = size(X); real = isreal(X); end
    tmp = randn(size_x)+(~real*randn(size_x)*1i);
    if isstruct(X), N = X; N.val = tmp;
    else N = tmp; end
    scale = frob(X)*10^(-SNR/10)/frob(N);
    if isstruct(X), X.val = X.val+scale*N.val;
    else X = X+scale*N; end
end

end
