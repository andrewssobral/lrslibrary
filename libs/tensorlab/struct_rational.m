function [x,state] = struct_rational(z,task,t)
%STRUCT_RATIONAL Matrix with columns as rational functions.
%   [x,state] = struct_rational(z,[],t) computes a matrix x in which the
%   jth column is equal to the rational function
%
%      polyval(z{1}(j,:),s)./polyval([z{2}(j,:) 1],s)
%
%   evaluated at the points s, defined as
%
%      (t-0.5*(min(t)+max(t)))/(0.5*(max(t)-min(t))).
%
%   The degree of the numerator and denomator are equal to size(z{1},2)-1
%   and size(z{2},2), respectively. The structure state stores information
%   which is reused in computing the right and left Jacobian-vector
%   products.
%
%   struct_rational(z,task,t) computes the right or left Jacobian-vector
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
%   See also struct_poly, struct_rbf.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 2, task = []; end
if nargin < 3
    error('struct_rational:t','Please supply evaluation points.');
end

% Center and scale evaluations points to fit in [-1,1].
mx = max(t);
mn = min(t);
t = (t-0.5*(mn+mx))/(0.5*(mx-mn));

if isempty(task) || (isempty(task.l) && isempty(task.r))
    x = polyval(z{1},t(:));
    if ~isempty(z{2})
        state.denom = polyval([z{2} ones(size(z{2},1),1)],t(:));
        x = x./state.denom;
        state.deriv = -x./state.denom;
    else
        state = [];
    end
elseif ~isempty(task.r)
    x = polyval(task.r{1},t(:));
    if ~isempty(z{2})
        x = x./task.denom+ ...
            task.deriv.*polyval([task.r{2} zeros(size(z{2},1),1)],t(:));
    end
    state = [];
elseif ~isempty(task.l)
    x = {zeros(size(z{1})),zeros(size(z{2}))};
    tmp = ones(length(t),1);
    if isempty(z{2})
        for i = size(z{1},2):-1:1
            x{1}(:,i) = sum(bsxfun(@times,conj(tmp),task.l),1);
            tmp = tmp.*t(:);
        end
    else
        for i = size(z{1},2):-1:1
            x{1}(:,i) = sum(conj(bsxfun(@rdivide,tmp,task.denom)).* ...
                task.l,1);
            tmp = tmp.*t(:);
        end
        tmp = t(:);
        for i = size(z{2},2):-1:1
            x{2}(:,i) = sum(conj(bsxfun(@times,tmp,task.deriv)).*task.l,1);
            tmp = tmp.*t(:);
        end
    end
    state = [];
end

end

function y = polyval(p,t)
    y = p(:,ones(length(t),1)).';
    for i = 2:size(p,2)
        y = bsxfun(@plus,bsxfun(@times,y,t),p(:,i).');
    end
end
