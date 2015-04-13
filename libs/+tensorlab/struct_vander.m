function [x,state] = struct_vander(z,task,deg)
%STRUCT_VANDER Vandermonde matrix.
%   [x,state] = struct_vander(z) uses the generator vector z to compute a
%   square Vandermonde matrix x as
%
%      [z(:).^0 z(:).^1 z(:).^2 ... z(:).^(length(z)-1)]
%
%   The structure state stores information which is reused in computing the
%   right and left Jacobian-vector products.
%
%   [x,state] = struct_vander(z,[],deg) generates x as the column vector 
%   z(:) raised to the powers deg(1):deg(2). The vector deg may contain
%   negative integers.
%
%   struct_vander(z,task,deg) computes the right or left Jacobian-vector
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
%   See also struct_hankel, struct_toeplitz.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 2, task = []; end
if nargin < 3 || isempty(deg), deg = length(z)-1; end
if isscalar(deg), deg = [0 deg]; end

if isempty(task) || (isempty(task.l) && isempty(task.r))
    tmp = vander(z(:),deg-1);
    x = bsxfun(@times,tmp,z(:));
    state.deriv = bsxfun(@times,tmp,deg(1):deg(2));
elseif ~isempty(task.r)
    x = bsxfun(@times,task.deriv,task.r(:));
    state = [];
elseif ~isempty(task.l)
    x = sum(conj(task.deriv).*task.l,2);
    state = [];
end

end

function V = vander(g,deg)

if deg(1) > 0
    V = ones(size(g,1),deg(2)+1);
    V(:,2) = g;
    for i = 3:deg(2)+1, V(:,i) = V(:,i-1).*g; end
    V = V(:,deg(1)+1:deg(2)+1);
elseif deg(2) < 0
    V = ones(size(g,1),-deg(1)+1);
    V(:,end-1) = 1./g;
    for i = -deg(1)-1:-1:1, V(:,i) = V(:,i+1)./g; end
    V = V(:,1:-deg(1)+deg(2)+1);
else
    V = ones(size(g,1),deg(2)-deg(1)+1);
    for i = -deg(1):-1:1, V(:,i) = V(:,i+1)./g; end
    V(:,-deg(1)+2,:) = g;
    for i = -deg(1)+3:deg(2)-deg(1)+1, V(:,i) = V(:,i-1).*g; end
end

end
