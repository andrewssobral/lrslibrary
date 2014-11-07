function [x,state] = struct_normalize(z,task)
%STRUCT_NORMALIZE Normalize columns to unit norm.
%   [x,state] = struct_normalize(z,[]) normalizes the columns of z so that
%   x(:,i) = z(:,i)/norm(z(:,i)) for all indices i. The structure state
%   stores information which is reused in computing the right and left
%   Jacobian-vector products.
%
%   struct_normalize(z,task) computes the right or left Jacobian-vector
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
%   See also struct_orth.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 2, task = []; end
state = [];

if isempty(task) || (isempty(task.l) && isempty(task.r))
    state.norm = dot(z,z,1);
    state.norm3 = state.norm.^(3/2);
    state.norm = sqrt(state.norm);
    x = bsxfun(@rdivide,z,state.norm);
elseif ~isempty(task.r)
    if ~isreal(z) || ~isreal(task.r)
        error('struct_normalize:nonanalytic',['Nonanalytic objective ' ...
            'functions are currently not supported in sdf_nls, please ' ...
            'use sdf_minf instead.']);
    end
    x = bsxfun(@rdivide,task.r,task.norm)- ...
        bsxfun(@times,z,sum(z.*task.r,1)./task.norm3);
elseif ~isempty(task.l)
    x = bsxfun(@rdivide,task.l,task.norm)- ...
        bsxfun(@times,z,sum(real(conj(z).*task.l),1)./task.norm3);
end

end
