function [x,state] = struct_sum(z,task,dim)
%STRUCT_SUM Sum of elements.
%   [x,state] = struct_sum(z,[],dim) computes x as sum(z,dim), or as
%   sum(z(:)) if dim is the empty matrix [] or not supplied. The structure
%   state stores information which is reused in computing the right and
%   left Jacobian-vector products.
%
%   struct_sum(z,task,dim) computes the right or left Jacobian-vector
%   product of this transformation, depending on the structure task. Use
%   the structure state and add the field 'r' of the same shape as z or the
%   field 'l' of the same shape as x to obtain the structure task for
%   computing the right and left Jacobian-vector products
%   
%      (dF(:)/dz(:).')*task.r(:) and
%      (dF(:)/dz(:).')'*task.l(:) + conj((dF(:)/dconj(z(:)).')'*task.l(:)),
%   
%   respectively. Here, F(z) represents this transormation, (:) signifies
%   vectorization and the derivative w.r.t. z (conj(z)) is a partial
%   derivative which treats conj(z) (z) as constant. The output has the
%   same shape as x or z for the right and left Jacobian-vector products,
%   respectively.
%   
%   See also struct_prod.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 2, task = []; end
if nargin < 3, dim = []; end
state = [];

if isempty(task)
    if isempty(dim), x = sum(z(:));
    else x = sum(z,dim);
    end
elseif ~isempty(task.r)
    if isempty(dim), x = sum(task.r(:));
    else x = sum(task.r,dim);
    end
elseif ~isempty(task.l)
    if isempty(dim), siz = size(z);
    else siz = ones(1,ndims(z)); siz(dim) = size(z,dim);
    end
    x = repmat(task.l,siz);
end

end
