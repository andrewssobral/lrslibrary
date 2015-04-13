function [x,state] = struct_band(z,task,size_mat,band)
%STRUCT_BAND Band matrix.
%   [x,state] = struct_band(z,[],size_mat,band) generates a band matrix x
%   of size size_mat, where the diagonals band(1) to band(2) are filled
%   column by column with entries from the vector z. For example, if x is a
%   square matrix of order n, the vector z should have length
%   sum(n-abs(band(1):band(2))). The structure state stores information
%   which is reused in computing the right and left Jacobian-vector
%   products.
%
%   struct_band(z,task,size_mat,band) computes the right or left
%   Jacobian-vector product of this transformation, depending on the
%   structure task. Use the structure state and add the field 'r' of the
%   same shape as z or the field 'l' of the same shape as x to obtain the
%   structure task for computing the right and left Jacobian-vector
%   products
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
%   See also struct_diag, struct_tridiag, struct_tril, struct_triu.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 2, task = []; end
if nargin < 3 || ~isvector(size_mat)
    error('struct_band:size_mat','Missing integer matrix order.');
end
if nargin < 4 || ~isvector(band)
    error('struct_band:band','Missing definition of band.');
end

if isempty(task) || (isempty(task.l) && isempty(task.r))
    state.idx = bsxfun(@plus,(size_mat(1)-1:-1:0).', ...
        (0:size_mat(2)-1)-size_mat(1)+1);
    state.idx = find(state.idx >= band(1) & state.idx <= band(2));
    x = zeros(size_mat);
    x(state.idx) = z;
elseif ~isempty(task.r)
    x = zeros(size_mat);
    x(task.idx) = task.r;
    state = [];
elseif ~isempty(task.l)
    x = task.l(task.idx);
    state = [];
end

end
