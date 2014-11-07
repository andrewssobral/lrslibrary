function [U,output] = cpd3_sgsd(T,U0,options)
%CPD3_SGSD CPD by simultaneous generalized Schur decomposition.
%   [U,output] = cpd3_sgsd(T,U0) computes the factor matrices U{1}, U{2}
%   and U{3} belonging to a canonical polyadic decomposition of the third-
%   order tensor T. The algorithm is initialized with the factor matrices
%   U0{n}. The number of rank-one terms R should satisfy sum(R <= size(T))
%   >= 2. The structure output returns additional information:
%
%      output.fval       - The value of the objective function, which
%                          is the Frobenius norm of the lower triangular
%                          entries of the 'upper triangular' matrices
%                          Rk(:,:,k) == Q*T(:,:,k)*Z, in every iteration.
%      output.info       - The circumstances under which the procedure
%                          terminated:
%                             1: Objective function tolerance reached.
%                             2: Maximum number of iterations reached.
%      output.iterations - The number of iterations.
%
%   cpd3_sgsd(T,U0,options) may be used to set the following options:
%
%      options.MaxIter = 200 - The maximum number of iterations.
%      options.TolFun = 1e-4 - The tolerance for the objective function
%                              between two successive iterates.

%   Authors: Mikael Sorensen (Mikael.Sorensen@kuleuven-kulak.be),
%            Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. De Lathauwer, B. De Moor, J. Vandewalle, "Computation of the
%       Canonical Decomposition by Means of a Simultaneous Generalized 
%       Schur Decomposition," SIAM J. Matrix Anal. Appl., Vol. 26, 2004,
%       pp. 295-327.
%   [2] A.-J. van der Veen, A. Paulraj, "An Analytical Constant Modulus 
%       Algorithm," IEEE Trans. Signal Proc., Vol. 44, 1996, pp. 1136-1155.

% Incomplete/sparse tensors not supported yet.
T = fmt(T,true);
if isstruct(T)
    error('This method does not yet support incomplete/sparse tensors.');
end

% Check the tensor T.
N = ndims(T);
if N ~= 3, error('cpd3_sgsd:T','ndims(T) should be 3.'); end

% Check the initial factor matrices U0.
U0 = U0(:).';
R = size(U0{1},2);
if any(cellfun('size',U0,2) ~= R)
    error('cpd3_sgsd:U0','size(U0{n},2) should be the same for all n.');
end
if sum(R > size(T)) > 1
    error('cpd3_sgsd:U0','sum(size(U0{n},2) <= size(T)) should be >= 2.');
end
if any(cellfun('size',U0,1) ~= size(T))
    error('cpd3_sgsd:U0','size(T,n) should equal size(U0{n},1).');
end

% Check the options structure.
if nargin < 3, options = struct; end
if ~isfield(options,'MaxIter'), options.MaxIter = 200; end
if ~isfield(options,'TolFun'), options.TolFun = 1e-4; end

% Prepermute the tensor so the longest two modes are first.
if size(T,1) < R || size(T,2) < R
    [~,idx] = sort(size(T),'descend');
    T = permute(T,idx);
    U0 = U0(idx);
    iperm(idx) = 1:N;
end

% Initialize Q and Z.
[Q,~] = qr(U0{1});
Q = Q';
[Z,~] = qr(fliplr(U0{2}));
Z = conj(fliplr(Z));

% Set up the algorithm.
[I,J,K] = size(T);
if I == R, Imax = R-1; else Imax = R; end
if J == R, Jmin = 2;   else Jmin = 1; end
T1 = reshape(T,I,[]);

% Simultaneous generalized Schur decomposition.
Rk = permute(reshape(reshape(permute(reshape(Q*T1,[I J K]), ...
     [1 3 2]),[],J)*Z,[I K J]),[1 3 2]);
output.fval = norm(tril(sum(abs(Rk).^2,3),J-R-1),'fro');
output.info = false;
output.iterations = 0;
while ~output.info

    % Update Q.
    Rk = reshape(Rk,I,[]);
    q = eye(I);
    for i = 1:Imax
        [u,~,~] = svd(q(i:end,:)*Rk(:,J-R+i:J:end));
        q(i:end,:) = u'*q(i:end,:);
    end
    Q = q*Q;
    
    % Update Z.
    Rk = q*Rk;
    z = eye(J);
    for j = R:-1:Jmin
        [~,~,v] = svd(reshape(Rk(j,:),J,K).'*z(:,1:J-R+j));
        z(:,1:J-R+j) = z(:,1:J-R+j)*v(:,end:-1:1);
    end
    Z = Z*z;
    
    % Apply Q and Z to Tk.
    Rk = permute(reshape(reshape(permute(reshape(Q*T1,[I J K]), ...
         [1 3 2]),[],J)*Z,[I K J]),[1 3 2]);
  
    % Update the output structure.
    output.fval(end+1) = norm(tril(sum(abs(Rk).^2,3),J-R-1),'fro');
    output.iterations = output.iterations+1;
    if abs(output.fval(end-1)-output.fval(end)) <= options.TolFun
        output.info = 1;
    end
    if output.iterations >= options.MaxIter, output.info = 2; end
    
end

% Reconstruct the factor matrices A, B and C.
Rk = Rk(1:R,end-R+1:end,:);
R1  = eye(R,R);
R11 = eye(R,R);
for i = R-1:-1:1
    for j = i+1:R
        A1 = Rk(j,j,:);
        A2 = Rk(i,i,:);
        b = Rk(i,j,:);
        for k = 1:K
            p = i+1:j-1;
            dk = diag(Rk(:,:,k));
            b(k) = b(k)-sum(R1(i,p).'.*dk(p).*R11(p,j));
        end
        rij = [A1(:) A2(:)]\b(:);
        R1(i,j)  = rij(1);
        R11(i,j) = rij(2);
    end
end
A = Q(1:R,:)'*R1;
B = (R11*Z(:,end-R+1:end)').';
C = reshape(T,[],K).'*conj(kr(B,A))/conj((B'*B).*(A'*A));
U = {A B C};

% Inverse permute the factor matrices.
if exist('iperm','var')
    U = U(iperm);
end
