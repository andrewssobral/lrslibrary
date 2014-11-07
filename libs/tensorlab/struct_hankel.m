function [x,state] = struct_hankel(z,task,size_mat,pre,post)
%STRUCT_HANKEL Hankel matrix.
%   [x,state] = struct_hankel(z) generates x as a square Hankel matrix in
%   which the elements of the vector z are placed on the anti-diagonals of
%   x from left to right. The structure state stores information which is
%   reused in computing the right and left Jacobian-vector products.
%
%   [x,state] = struct_hankel(z,[],size_mat,pre,post) sets the size of the
%   Hankel matrix equal to size_mat and uses the elements of the vector
%   [pre(:);z(:);post(:)] for the anti-diagonals of x.
%
%   struct_hankel(z,task,size_mat,pre,post) computes the right or left
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
%   See also struct_toeplitz, struct_vander.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 2, task = []; end
if nargin < 3 || isempty(size_mat), size_mat = (length(z)+1)/2*[1 1]; end
if nargin < 4, pre = []; end
if nargin < 5, post = []; end
state = [];

if isempty(task) || (isempty(task.l) && isempty(task.r))
    g = [pre(:);z(:);post(:)];
    x = hankel(g(1:size_mat(1)),g(size_mat(1):end));
elseif ~isempty(task.r)
    g = [zeros(numel(pre),1);task.r(:);zeros(numel(post),1)];
    x = hankel(g(1:size_mat(1)),g(size_mat(1):end));
elseif ~isempty(task.l)
    task.l = fliplr(task.l);
    x = arrayfun(@(d)sum(diag(task.l,d)), ...
        (size_mat(2)-1-numel(pre):-1:-size_mat(1)+1+numel(post))');
end

end
