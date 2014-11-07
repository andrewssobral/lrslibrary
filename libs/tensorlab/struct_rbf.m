function [x,state] = struct_rbf(z,task,t)
%STRUCT_RBF Matrix with columns as sums of Gaussian RBF kernels.
%   [x,state] = struct_rbf(z,[],t) computes a matrix x in which the
%   jth column is equal to
%
%      sum_{n=1}^N z{1}(j,n)*exp(-(t-z{2}(j,n)).^2./(2*z{3}(j,n)^2))
%
%   evaluated at the points t. The size of each cell in z should be equal
%   to size(x,2)-by-N. The structure state stores information which is
%   reused in computing the right and left Jacobian-vector products.
%
%   struct_rbf(z,task,t) computes the right or left Jacobian-vector product
%   of this transformation, depending on the structure task. Use the
%   structure state and add the field 'r' of the same shape as z or the 
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
%   See also struct_poly, struct_rational.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 2, task = []; end
if nargin < 3
    error('struct_rbf:t','Please supply evaluation points.');
end

if isempty(task) || (isempty(task.l) && isempty(task.r))
    tmp = bsxfun(@minus,t(:),permute(z{2},[3 1 2]));
    state.ic2 = permute(1./(z{3}.^2),[3 1 2]);
    state.ic3 = state.ic2./permute(z{3},[3 1 2]);
    state.exp = exp(-bsxfun(@times,tmp.^2,0.5*state.ic2));
    state.g = bsxfun(@times,permute(z{1},[3 1 2]),state.exp);
    state.gt = state.g.*tmp;
    state.gt2 = state.gt.*tmp;
    state.gtic2 = bsxfun(@times,state.gt,state.ic2);
    state.gt2ic3 = bsxfun(@times,state.gt2,state.ic3);
    x = sum(state.g,3);
elseif ~isempty(task.r)
    x = sum(bsxfun(@times,task.exp,permute(task.r{1},[3 1 2]))+ ...
        bsxfun(@times,task.gt,task.ic2.*permute(task.r{2},[3 1 2]))+ ...
        bsxfun(@times,task.gt2,task.ic3.*permute(task.r{3},[3 1 2])),3);
    state = [];
elseif ~isempty(task.l)
    x = cell(1,3);
    x{1} = permute(sum(bsxfun(@times,conj(task.exp),task.l),1),[2 3 1]);
    x{2} = permute(sum(bsxfun(@times,conj(task.gtic2),task.l),1),[2 3 1]);
    x{3} = permute(sum(bsxfun(@times,conj(task.gt2ic3),task.l),1),[2 3 1]);
    state = [];
end

end
