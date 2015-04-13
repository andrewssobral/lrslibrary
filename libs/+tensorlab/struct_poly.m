function [x,state] = struct_poly(z,task,t)
%STRUCT_POLY Matrix with columns as polynomials.
%   [x,state] = struct_poly(z,[],t) computes a matrix x in which the
%   jth column is equal to the polynomial
%
%      polyval(z(j,:),s)
%
%   evaluated at the points s, defined as
%
%      (t-0.5*(min(t)+max(t)))/(0.5*(max(t)-min(t))).
%
%   The degree of the polynomial is equal to size(z,2)-1. The structure
%   state stores information which is reused in computing the right and
%   left Jacobian-vector products.
%
%   struct_poly(z,task,t) computes the right or left Jacobian-vector
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
%   See also struct_rational, struct_rbf.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 2, task = []; end
if nargin < 3
    error('struct_poly:t','Please supply evaluation points.');
end

if isempty(task)
    [x,state] = struct_rational({z,[]},task,t);
elseif ~isempty(task.r)
    task.r = {task.r,[]};
    [x,state] = struct_rational({z,[]},task,t);
elseif ~isempty(task.l)
    [x,state] = struct_rational({z,[]},task,t); x = x{1};
end

end
