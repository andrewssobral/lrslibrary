function [x,state] = struct_matvec(z,task,mat,mat2)
%STRUCT_MATVEC Matrix-vector and matrix-matrix product.
%   [x,state] = struct_matvec(z,[],mat) computes x as the matrix-vector (or
%   matrix-matrix) product mat*z. The structure state stores information
%   which is reused in computing the right and left Jacobian-vector
%   products.
%
%   [x,state] = struct_matvec(z,[],mat,mat2) computes x as the
%   matrix-matrix product mat*z*mat2. Alternatively, set mat equal to the
%   empty matrix [] to compute x as z*mat2.
%
%   struct_matvec(z,task,mat,mat2) computes the right or left
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

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. Sorber, M. Van Barel, L. De Lathauwer, "Structured data fusion,"
%       ESAT-SISTA Internal Report 13-177, KU Leuven, 2013.

if nargin < 2, task = []; end
if nargin < 4, mat2 = []; end
state = [];

if isempty(task) || (isempty(task.l) && isempty(task.r))
    if isempty(mat2)
        x = mat*z;
    elseif isempty(mat)
        x = z*mat2;
    else
        x = mat*z*mat2;
    end
elseif ~isempty(task.r)
    if isempty(mat2)
        x = mat*task.r;
    elseif isempty(mat)
        x = task.r*mat2;
    else
        x = mat*task.r*mat2;
    end
elseif ~isempty(task.l)
    if isempty(mat2)
        x = mat'*task.l;
    elseif isempty(mat)
        x = task.l*mat2';
    else
        x = mat'*task.l*mat2';
    end
end

end
