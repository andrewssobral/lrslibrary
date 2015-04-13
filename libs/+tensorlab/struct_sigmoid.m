function [x,state] = struct_sigmoid(z,task,range)
%STRUCT_SIGMOID Constrain array elements to an interval.
%   [x,state] = struct_sigmoid(z,[],range) computes x by applying a scaled
%   and translated sigmoid transformation to the array z such that
%   all(x(:) > range(1)) && all(x(:) < range(2)). The structure state
%   stores information which is reused in computing the right and left
%   Jacobian-vector products.
%
%   struct_sigmoid(z,task,range) computes the right or left Jacobian-vector
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
%   See also struct_abs, struct_nonneg.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 2, task = []; end
if nargin < 3, range = [-1 1]; end

if isempty(task) || (isempty(task.l) && isempty(task.r))
    tmp = 1+z.^2;
    tmp2 = (0.5*(range(2)-range(1)));
    x = tmp2*(z./sqrt(tmp)+1)+range(1);
    state.deriv = tmp2./(tmp.*sqrt(tmp));
elseif ~isempty(task.r)
    x = task.deriv.*task.r;
    state = [];
elseif ~isempty(task.l)
    x = conj(task.deriv).*task.l;
    state = [];
end

end
