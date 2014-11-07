function E = lmlrares(T,U,S,options)
%LMLRARES Residual of a LMLRA.
%   E = lmlrares(T,U,S) computes the residual tensor E as lmlragen(U,S)-T.
%
%   See also btdres, cpdres.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check the options structure.
if nargin < 4, options = struct; end

% Compute the residual.
E = btdres(T,{[U(:).',S]},options);
