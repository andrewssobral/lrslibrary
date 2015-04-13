function [x,state] = struct_power(z,task,deg)
%STRUCT_POWER Array power.
%   [x,state] = struct_power(z,[],deg) computes x as z.^deg. The structure
%   state stores information which is reused in computing the right and
%   left Jacobian-vector products.
%
%   struct_power(z,task,deg) computes the right or left Jacobian-vector
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
%   See also struct_nonneg, struct_sigmoid.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 2, task = []; end
if nargin < 3, deg = 2; end

if isempty(task) || (isempty(task.l) && isempty(task.r))
    if isscalar(deg) && deg == 1, tmp = 1;
    elseif isscalar(deg) && deg == 2, tmp = z;
    else tmp = z.^(deg-1);
    end
    x = z.*tmp;
    state.deriv = deg.*tmp;
elseif ~isempty(task.r)
    x = task.deriv.*task.r;
    state = [];
elseif ~isempty(task.l)
    x = conj(task.deriv).*task.l;
    state = [];
end

end
