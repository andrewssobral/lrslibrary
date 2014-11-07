function E = cpdres(T,U,options)
%CPDRES Residual of a polyadic decomposition.
%   E = cpdres(T,U) computes the residual tensor E as cpdgen(U)-T.
%
%   See also btdres, lmlrares.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check the options structure.
if nargin < 3, options = struct; end
if ~isfield(options,'TolLargeScale'), options.TolLargeScale = 0.02; end

% Compute the residual.
T = fmt(T,true);
if isstruct(T)
    if ~T.incomplete || length(T.ind)/prod(T.size) > options.TolLargeScale
        E = U{1}*kr(U(end:-1:2)).';
        if T.incomplete, E = E(T.ind)-T.val;
        elseif T.sparse, E(T.ind) = E(T.ind)-T.val;
        else E = E-reshape(T,size(E));
        end
    else
        E = U{1}(T.sub{1},:);
        for n = 2:length(U)
            E = E.*U{n}(T.sub{n},:);
        end
        E = sum(E,2)-T.val;
    end
else
    E = cpdgen(U)-T;
end
