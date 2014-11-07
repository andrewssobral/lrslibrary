function [x,state] = struct_tridiag(z,task)
%STRUCT_TRIDIAG Tridiagonal matrix.
%   [x,state] = struct_tridiag(z) generates x as a tridiagonal matrix by
%   using the vector z to fill the matrix column by column. For a matrix x
%   of order n, the vector z should have length 3*n-2. The structure state
%   stores information which is reused in computing the right and left
%   Jacobian-vector products.
%
%   struct_tridiag(z,task) computes the right or left Jacobian-vector
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
%   See also struct_band, struct_diag, struct_tril, struct_triu.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 2, task = []; end
n = (length(z)+2)/3;
[x,state] = struct_band(z,task,[n n],[-1 1]);

end
