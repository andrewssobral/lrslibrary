function [x,state] = struct_orth(z,task,size_mat)
%STRUCT_ORTH Rectangular matrix with orthonormal columns.
%   [x,state] = struct_orth(z,[],size_mat) computes x as a rectangular
%   matrix of size size_mat with orthonormal columns, given the vector z.
%   The vector z should have length size_mat(2)*(size_mat(1)- ...
%   (size_mat(2)-1)/2). The structure state stores information which is
%   reused in computing the right and left Jacobian-vector products.
%
%   struct_orth(z,task,size_mat) computes the right or left Jacobian-vector
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
%   See also struct_normalize.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 2, task = []; end
N = length(z);
if nargin < 3, size_mat = 0.5*(-1+sqrt(1+8*N))*[1 1]; end
m = size_mat(1);
n = size_mat(2);

if isempty(task) || (isempty(task.l) && isempty(task.r))
    
    % Initialize state.
    state.offset = cumsum([1 (m-n+1):m]);
    state.left = cell(1,n);
    state.right = cell(1,n);
    state.projleft = cell(1,n);
    state.projright = cell(1,n);
    state.ivv = 1./arrayfun(@(n)z(state.offset(n):state.offset(n+1)-1)'*...
        z(state.offset(n):state.offset(n+1)-1),1:n);
    
    % Apply product of Householder reflectors from right to left.
    % Used for right Jacobian-vector product.
    v = z(state.offset(1):state.offset(2)-1);
    state.right{1} = eye(length(v),1)-v*(conj(v(1))*2*state.ivv(1));
    for i = 2:n
        v = z(state.offset(i):state.offset(i+1)-1);
        tmpa = -v*(conj(v(1))*2*state.ivv(i)); tmpa(1) = tmpa(1)+1;
        tmpb = state.right{i-1}; tmpc = [zeros(1,size(tmpb,2));tmpb];
        state.right{i} = [tmpa tmpc-((2*state.ivv(i))*v)*(v(2:end)'*tmpb)];
    end
    x = state.right{end};
    
    % Apply product of Householder reflectors from left to right.
    % Used for right Jacobian-vector product.
    if nargout > 1
    v = z(state.offset(end-1):state.offset(end)-1);
    state.left{n} = eye(m)-v*(v'*(2*state.ivv(n)));
    for i = n-1:-1:2
        v = z(state.offset(i):state.offset(i+1)-1);
        tmp = state.left{i+1}(:,end-length(v)+1:end);
        state.left{i} = tmp-(tmp*v)*(v'*(2*state.ivv(i)));
        state.left{i+1} = state.left{i+1}(:,2:end);
    end
    if n > 1, state.left{2} = state.left{2}(:,2:end); end
    
    % Cache projected vectors. Used for left Jacobian-vector product.
    for i = 1:n
        v = z(state.offset(i):state.offset(i+1)-1);
        if i < n
            state.projleft{i} = state.left{i+1}*v;
        else
            state.projleft{i} = v;
        end
        if i > 1
            state.projright{i} = [conj(v(1)) v(2:end)'*state.right{i-1}];
        else
            state.projright{i} = conj(v(1));
        end
    end
    end

elseif ~isempty(task.r)
    
    % Compute J(z)*r as a sum of matrices of the form
    % [H1*...*Hn-1*(dHn/dvn*rn)*Hn+1*...*Hn].
    if ~isreal(z) || ~isreal(task.r)
        error('struct_orth:nonanalytic',['Nonanalytic objective ' ...
            'functions are currently not supported in sdf_nls, please ' ...
            'use sdf_minf instead.']);
    end
    v = z(task.offset(1):task.offset(2)-1); ivv = task.ivv(1);
    r = task.r(task.offset(1):task.offset(2)-1);
    x = r*(v(1)*(-2*ivv))+v*(r(1)*(-2*ivv)+v(1)*(4*ivv*ivv*sum(v.*r)));
    if n >= 2
        x = [zeros(m,m-size(x,1)) task.left{2}*x];
        for i = 2:n
            v = z(task.offset(i):task.offset(i+1)-1); ivv = task.ivv(i);
            r = task.r(task.offset(i):task.offset(i+1)-1);
            tmp = [v'*(-2*ivv); r'*(-2*ivv)+v'*(4*ivv*ivv*sum(v.*r))];
            H = [tmp(:,1) tmp(:,2:end)*task.right{i-1}];
            if i < n, H = (task.left{i+1}*[r v])*H;
            else H = [r v]*H;
            end
            x(:,end-size(H,2)+1:end) = x(:,end-size(H,2)+1:end)+H;
        end
    end
    state = [];
    
elseif ~isempty(task.l)
    
    % Compute (dF/dz^T)'*l including a small part of conj((dF/dz^H)'*l).
    x = zeros(size(z));
    for i = 1:n
        tmp = task.l(:,end-i+1:end)*task.projright{i}';
        projvv = task.projleft{i}'*tmp;
        if i < n
            projeye = task.left{i+1}'*tmp;
        else
            projeye = tmp;
        end
        idx = task.offset(i):task.offset(i+1)-1;
        x(idx) = x(idx)-(2*task.ivv(i))*projeye+ ...
            2*real(2*task.ivv(i)^2*projvv).*z(idx);
    end
    
    % Compute the remainder of conj((dF/dz^H)'*l).
    for i = 1:n
        tmp = task.projleft{i}'*task.l;
        if i > 1
            projeye = [tmp(end-i+1) tmp(:,end-i+2:end)*task.right{i-1}'];
        else
            projeye = [tmp(end) zeros(1,task.offset(2)-2)];
        end
        idx = task.offset(i):task.offset(i+1)-1;
        x(idx) = x(idx)-(2*task.ivv(i))*projeye';
    end
    state = [];
    
end

end
