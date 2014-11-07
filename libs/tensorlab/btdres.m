function E = btdres(T,U,options)
%BTDRES Residual of a BTD.
%   E = btdres(T,U) computes the residual tensor E as btdgen(U)-T.
%
%   See also cpdres, lmlrares.

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
        E = U{1}{1}*mtkronprod(U{1}{end},U{1}(1:end-1),1,'H');
        for r = 2:length(U)
            E = E+U{r}{1}*mtkronprod(U{r}{end},U{r}(1:end-1),1,'H');
        end
        if T.incomplete, E = E(T.ind)-T.val;
        elseif T.sparse, E(T.ind) = E(T.ind)-T.val;
        else E = E-reshape(T,size(E));
        end
    else
        E = -T.val;
        for r = 1:length(U)
            S = U{r}{end};
            size_core = cellfun('size',U{r}(1:end-1),2);
            idx = cell(1,length(U{r})-1);
            for i = 1:numel(S)
                [idx{:}] = ind2sub(size_core,i);
                tmp = S(idx{:})*U{r}{1}(T.sub{1},idx{1});
                for n = 2:length(size_core)
                    tmp = tmp.*U{r}{n}(T.sub{n},idx{n});
                end
                E = E+tmp;
            end
        end
    end
else
    E = btdgen(U)-T;
end
