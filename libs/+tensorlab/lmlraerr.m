function sangle = lmlraerr(U,Uest)
%LMLRAERR Errors between factor matrices in a LMLRA.
%   sangle = lmlraerr(U,Uest), with U and Uest cell arrays of factor
%   matrices, computes the angle sangle(n) between the subspaces U{n} and
%   Uest{n} for n = 1:length(U). The factor matrices U{n} and Uest{n}
%   should have an equal number of rows, but the number of columns may
%   differ.
%
%   See also lmlragen, cpderr.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check input.
if ~iscell(U), U = {U}; end;
if ~iscell(Uest), Uest = {Uest}; end;    
if length(U) ~= length(Uest)
    error('lmlraerr:U','length(U) should equal length(Uest).');
end
if any(cellfun('size',U(:),1) ~= cellfun('size',Uest(:),1))
    error('lmlraerr:U','size(U{n},1) should equal size(Uest{n},1).');
end

% Compute angle between mode-n subspaces.
sangle = cellfun(@(u,v)subspace(u,v),U(:).',Uest(:).');
