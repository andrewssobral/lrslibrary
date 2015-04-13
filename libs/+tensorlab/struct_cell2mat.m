function [x,state] = struct_cell2mat(z,task)
%STRUCT_CELL2MAT Convert the contents of a cell array into a matrix.
%   [x,state] = struct_cell2mat(z) generates the matrix x as the function
%   cell2mat applied to the cell array z. The structure state stores
%   information which is reused in computing the right and left
%   Jacobian-vector products.
%
%   struct_cell2mat(z,task) computes the right or left Jacobian-vector
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

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 2, task = []; end

if isempty(task) || (isempty(task.l) && isempty(task.r))
    tmp = cellfun(@size,z,'UniformOutput',false);
    N = max(cellfun(@length,tmp));
    state.dim = cell(1,N);
    for n = 1:N
        s = size(tmp,n);
        idx = mat2cell([ones(s,n-1) (1:s)' ones(s,N-n)],s,ones(1,N));
        idx = sub2ind(size(tmp),idx{:});
        state.dim{n} = cellfun(@(s)s(n),tmp(idx));
    end
    x = cell2mat(z);
elseif ~isempty(task.r)
    x = cell2mat(task.r);
    state = [];
elseif ~isempty(task.l)
    x = mat2cell(task.l,task.dim{:});
    state = [];
end

end
