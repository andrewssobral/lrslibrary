function [U,S,output] = lmlra3_dgn(T,U0,options)
%LMLRA3_DGN LMLRA by a differential-geometric Newton method.
%   [U,S,output] = lmlra3_dgn(T,U0) computes the factor matrices U{1}, 
%   U{2} and U{3} and core tensor S belonging to a low multilinear rank
%   approximation of the third-order tensor T. The algorithm is initialized
%   with the factor matrices U0{n}. Global convergence is not guaranteed,
%   hence it is recommended to choose U0 close to a local minimum, e.g. by
%   initializing with a truncated MLSVD. The structure output returns
%   additional information:
%
%      output.info       - The circumstances under which the procedure
%                          terminated:
%                             1: Objective function tolerance reached.
%                             2: Maximum number of iterations reached.
%      output.iterations - The number of iterations.
%      output.normF      - The norm of the residuals F in every iteration.
%      output.sangle     - The difference in subspace angle of U{1} between
%                          every two successive iterates.
%
%   lmlra3_dgn(T,U0,options) may be used to set the following options:
%
%      options.BiCGTol = 1e-6 - The tolerance for the Newton equation
%                               solved by BiCGSTAB.
%      options.MaxIter = 500  - The maximum number of iterations.
%      options.TolFun = 1e-4  - The tolerance for the norm of the
%                               residuals F, where F is defined as in [1].

%   Authors: Mariya Ishteva (mariya.ishteva@cc.gatech.edu)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] M. Ishteva, L. De Lathauwer, P.-A. Absil, S. Van Huffel,
%       "Differential-geometric Newton method for the best rank-(R1,R2,R3)
%       approximation of tensors", Numerical Algorithms, Vol. 51, No. 2,
%       2009, pp. 179-194.

% Check the tensor T.
N = ndims(T);
if N ~= 3, error('lmlra3_dgn:T','ndims(T) should be 3.'); end
if ~isreal(T), error('lmlra3_dgn:T','isreal(T) should be true.'); end

% Check the initial factor matrices U0.
if any(cellfun('size',U0(:).',1) ~= size(T))
    error('lmlra3_dgn:U0','size(T,n) should equal size(U0{n},1).');
end

% Check the options structure.
if nargin < 3, options = struct; end
if ~isfield(options,'MaxIter'), options.MaxIter = 100; end
if ~isfield(options,'TolFun'), options.TolFun = 1e-4; end
if ~isfield(options,'BiCGTol'), options.BiCGTol = 1e-6; end

% Initialization
[I1,R1] = size(U0{1});
[I2,R2] = size(U0{2});
[I3,R3] = size(U0{3});
Xnew = U0;
Z = {zeros(I1,R1),zeros(I2,R2),zeros(I3,R3)};
%vecZ = Z(:);
niter = 0;
AX = struct;

% Dimension of the quotient manifold R^{nxp}_*/O(p):
dim1 = I1*R1 - R1*(R1-1)/2;
dim2 = I2*R2 - R2*(R2-1)/2;
dim3 = I3*R3 - R3*(R3-1)/2;

% Cost function at U0 (and some additional variables)
[F,AX] = tensor_F(T,U0,AX);
normF = norm(F{1}(:)) + norm(F{2}(:)) + norm(F{3}(:));

%warning('off','MATLAB:gmres:tooSmallTolerance');
output.info = false;
output.sangle = inf;
while ~output.info

    niter = niter+1;
    X = Xnew;

    % Additional variables
    % Note: The variables in the structure AX are computed once at the
    % beginning of the GMRES step to avoid multiple computations of the
    % expressions at each step within GMRES.
    AX.invUTU = (X{1}'*X{1})\eye(R1);
    AX.invVTV = (X{2}'*X{2})\eye(R2);
    AX.invWTW = (X{3}'*X{3})\eye(R3);
    AX.invUTU_RURU = AX.invUTU*AX.RURU;
    AX.invVTV_RVRV = AX.invVTV*AX.RVRV;
    AX.invWTW_RWRW = AX.invWTW*AX.RWRW;
    AX.skew_invUTU_RURU = skew(AX.invUTU_RURU);
    AX.skew_invVTV_RVRV = skew(AX.invVTV_RVRV);
    AX.skew_invWTW_RWRW = skew(AX.invWTW_RWRW);

    % The operator (LHS of the system)
    M_op = @(vecZ)tensor_vec_PhDPhF(T,X,vecZ,AX);

    % The operator, viewed as an operator from the horizontal
    % space onto the horizontal space, is nonsingular in general.
    % Hence, we use GMRES to solve the Newton equation. All the inner
    % iterates produced by GMRES are in the horizontal space, since
    % the linear operator is onto the horizontal space.

    % The RHS of the system
    PhF = tensor_Ph(X,F);
    vecPhF = [PhF{1}(:); PhF{2}(:); PhF{3}(:)];

    % Main computation: solving the system using GMRES
    %[vecZ,flagGMRES] = gmres(M_op,-vecPhF,[],options.GMRESTol, ...
    %    dim1+dim2+dim3);
    %if flagGMRES ~= 0 && flagGMRES ~= 3
    %    warning('lmlra3_dgn:gmres','flag = %i',flagGMRES);
    %end
    [vecZ,flagBCG] = bicgstab(M_op,-vecPhF,options.BiCGTol,dim1+dim2+dim3);
    if flagBCG > 1, warning('lmlra3_dgn:bicgstab','flag = %i',flagBCG); end

    % Reshape the output vecZ:
    Z{1} = reshape(vecZ(1:I1*R1),I1,R1);
    Z{2} = reshape(vecZ(I1*R1+1:I1*R1+I2*R2),I2,R2);
    Z{3} = reshape(vecZ(I1*R1+I2*R2+1:end),I3,R3);

    % Apply the Newton update:
    Xnew{1} = X{1}+Z{1};
    Xnew{2} = X{2}+Z{2};
    Xnew{3} = X{3}+Z{3};

    % Cost function at Xnew (and some additional variables)
    [F,AX] = tensor_F(T,Xnew,AX);
    normF(end+1) = norm(F{1}(:))+norm(F{2}(:))+norm(F{3}(:));
    
    % Update the output structure.
    output.info = normF(end) <= options.TolFun;
    output.sangle(end+1) = subspace(Xnew{1},X{1});
    if niter >= options.MaxIter, output.info = 2; end
    
end

U = Xnew{1};
V = Xnew{2};
W = Xnew{3};

% Update the output structure.
output.iterations = niter;
output.normF = normF;

% Compute the final core tensor.
U = {U,V,W};
S = tmprod(T,U,[1 2 3],'H');

end

% *************************************
function S = skew(T)
S = .5*(T - T');
end

function [A1,A2,A3] = tens2mats(A)
% tens2mats matricizes a tensor.
%   tens2mats(A) matricizes a third-order tensor A w.r.t all three modes,
%   such that the columns of the output matrix are the mode-n fibers of the
%   tensor. The mode-n fibers are placed in the matrix according to the
%   ordering in
%
%   L. De Lathauwer, B. De Moor, J. Vandewalle, "A Multilinear Singular 
%   Value Decomposition", SIAM J. Matrix Anal. Appl., Vol. 21, No. 4, 2000,
%   pp. 1253-1278.

% Version: 30.08.2010

% called functions: --
 
[I,J,K] = size(A);

A1 = reshape(permute(A,[1 3 2]),[I J*K]);
A2 = reshape(permute(A,[2 1 3]),[J K*I]);
A3 = reshape(permute(A,[3 2 1]),[K I*J]);
end

function [F,AX]  = tensor_F(A,X,AX)
% tensor_F(A,X,AX) computes the cost function F of the geometric Newton
% algorithm in
%
% M. Ishteva, L. De Lathauwer, P.-A. Absil, S. Van Huffel,
% "Differential-geometric Newton method for the best rank-(R1,R2,R3)
% approximation of tensors", Numerical Algorithms, Vol. 51, No. 2, 2009,
% pp. 179–194.
%
% (It also stores some of the intermediate results in the AX structure.)
%
% F = {F1 F2 F3}
% F1 = U*RU*RU' - A1*kron(V,W)*RU'
% F2 = V*RV*RV' - A2*kron(W,U)*RV'
% F3 = W*RW*RW' - A3*kron(U,W)*RW'
% where
% X = {U V W}, A1,A2,A3 are matricized versions of A, see tens2mats(A)
% RU = U'*A1*kron(V,W)
% RV = V'*A2*kron(W,U)
% RW = W'*A3*kron(U,V)

% Version: 30.08.2010

% called functions: tmprod, tens2mats

U = X{1};
V = X{2};
W = X{3};

[I1,R1] = size(U);
[I2,R2] = size(V);
[I3,R3] = size(W);

AX.AU = tmprod(A,U',1);
AX.AV = tmprod(A,V',2);
AX.AW = tmprod(A,W',3);

AX.AUV = tmprod(AX.AU,V',2);
AX.AVW = tmprod(AX.AV,W',3);
AX.AWU = tmprod(AX.AW,U',1);

AX.AUVW = tmprod(AX.AUV,W',3);

AX.A1VW = reshape(permute(AX.AVW,[1 3 2]),[I1 R2*R3]);
AX.A2WU = reshape(permute(AX.AWU,[2,1,3]),[I2 R3*R1]);
AX.A3UV = reshape(permute(AX.AUV,[3,2,1]),[I3 R1*R2]);
AX.A1VW_A1VW = AX.A1VW*AX.A1VW';
AX.A2WU_A2WU = AX.A2WU*AX.A2WU';
AX.A3UV_A3UV = AX.A3UV*AX.A3UV';

[AX.RU,AX.RV,AX.RW] = tens2mats(AX.AUVW);

AX.RURU = AX.RU*AX.RU';
AX.RVRV = AX.RV*AX.RV';
AX.RWRW = AX.RW*AX.RW';
AX.A1VW_RU = AX.A1VW*AX.RU';
AX.A2WU_RV = AX.A2WU*AX.RV';
AX.A3UV_RW = AX.A3UV*AX.RW';

F{1} = U*AX.RURU - AX.A1VW_RU;
F{2} = V*AX.RVRV - AX.A2WU_RV;
F{3} = W*AX.RWRW - AX.A3UV_RW;
end

function res = tensor_vec_PhDPhF(A,X,vecZ,AX)
% tensor_vec_PhDPhF(A,X,vecZ,AX) returns the vec of the horizontal
% projection of the directional derivative of the horizontal projection of
% F (see tensor_F.m) along the direction Z.
%
% The direction Z is meant to be horizontal in the sense of tensor_Ph, but
% this condition is not necessary for the call to succeed.
%
% This function is implemented for the purpose of being used in a
% call to solve the Newton equation Ph D (PhF)(X)[Z] = -PhF(X).

% Version: 30.08.2010

% called functions: tensor_Ph, tensor_DPhF

[I1,R1] = size(X{1});
[I2,R2] = size(X{2});
[I3,R3] = size(X{3});

Z1 = reshape(vecZ(1:I1*R1),I1,R1);
Z2 = reshape(vecZ(I1*R1+1:I1*R1+I2*R2),I2,R2);
Z3 = reshape(vecZ(I1*R1+I2*R2+1:end),I3,R3);
Z = {Z1,Z2,Z3};

PhDPhFXZ = tensor_Ph(X, tensor_DPhF(A,X,Z,AX));

res = [PhDPhFXZ{1}(:); PhDPhFXZ{2}(:); PhDPhFXZ{3}(:)];
end

function res = tensor_Ph(X,Z)
% tensor_Ph(X,Z) returns the projection onto the horizontal space
% H_X = {XS+X_perpK: S'=S}
% along the vertical space
% V_X = {XOmega: Omega'=-Omega}
% of the quotient manifold R*nxp/O(p).
% Ph(X,Z) = Z - X skew(inv(X'X)X'Z),
% where skew(B) = (B-B')/2.
%
% Note that the vertical and horizontal spaces are not orthogonal with
% respect to the canonical metric on R*nxp, but rather with respect to
% a certain modified metric.

% Version: 30.08.2010

% called functions: --

U = X{1};
V = X{2};
W = X{3};

S1 = (U'*U)\(U'*Z{1});
res{1} = Z{1}-U*(S1-S1')/2;

S2 = (V'*V)\(V'*Z{2});
res{2} = Z{2}-V*(S2-S2')/2;

S3 = (W'*W)\(W'*Z{3});
res{3} = Z{3}-W*(S3-S3')/2;
end

function res = tensor_DPhF(~,X,Z,AX)
% tensor_DPhF(A,X,Z,AX) returns the directional derivative of the
% horizontal projection of F (see tensor_F.m) along the direction Z.
%
% The direction Z is meant to be horizontal in the sense of tensor_Ph, but
% this condition is not necessary for the call to succeed.

% This is equation (12) in:
%
% M. Ishteva, L. De Lathauwer, P.-A. Absil, S. Van Huffel,
% "Differential-geometric Newton method for the best rank-(R1,R2,R3)
% approximation of tensors", Numerical Algorithms, 51(2):179–194, 2009.
%
% The computation is optimized so that no Kronecker products are computed.
% Part of the help variables are stored in AX.

% Version: 30.08.2010

% called functions: tmprod

U = X{1};
V = X{2};
W = X{3};

Z1 = Z{1};
Z2 = Z{2};
Z3 = Z{3};

[I1,R1] = size(U);
[I2,R2] = size(V);
[I3,R3] = size(W);

% 1
Z1A1VWRU = Z1' * AX.A1VW_RU ;

AZ2W = tmprod(AX.AW,Z2',2);
A1Z2W = reshape(permute(AZ2W,[1 3 2]),[I1 R2*R3]);
AVZ3 = tmprod(AX.AV,Z3',3);
A1VZ3 = reshape(permute(AVZ3,[1 3 2]),[I1 R2*R3]);

ADVDWAU = (A1Z2W + A1VZ3)*AX.RU' + AX.A1VW*(A1Z2W' + A1VZ3')*U;

arg_of_skew = AX.invUTU*( ...
    -(Z1'*U + U'*Z1)*AX.invUTU_RURU + Z1A1VWRU + Z1A1VWRU'+ U'*ADVDWAU);

res{1} = Z1*AX.RURU + U*(Z1A1VWRU + Z1A1VWRU') ...
    + U*(U'*ADVDWAU) - AX.A1VW_A1VW*Z1 - ADVDWAU ...
    + Z1*AX.skew_invUTU_RURU + U*skew(arg_of_skew);

% 2
Z2A2WURV = Z2'*AX.A2WU_RV;

AZ3U = tmprod(AX.AU,Z3',3);
A2Z3U = reshape(permute(AZ3U,[2,1,3]),[I2 R3*R1]);
AWZ1 = tmprod(AX.AW,Z1',1);
A2WZ1 = reshape(permute(AWZ1,[2,1,3]),[I2 R3*R1]);

ADWDUAV = (A2Z3U + A2WZ1)*AX.RV' + AX.A2WU*(A2Z3U' + A2WZ1')*V;

arg_of_skew = AX.invVTV*( ...
    -(Z2'*V + V'*Z2)*AX.invVTV_RVRV + Z2A2WURV + Z2A2WURV' + V'*ADWDUAV);

res{2} = Z2*AX.RVRV + V*(Z2A2WURV + Z2A2WURV') ...
    + V*(V'*ADWDUAV) - AX.A2WU_A2WU*Z2 - ADWDUAV ...
    + Z2*AX.skew_invVTV_RVRV + V*skew(arg_of_skew);

% 3
Z3A3UVRW = Z3'*AX.A3UV_RW;

AZ1V = tmprod(AX.AV,Z1',1);
A3Z1V = reshape(permute(AZ1V,[3,2,1]),[I3 R1*R2]);
AUZ2 = tmprod(AX.AU,Z2',2);
A3UZ2 = reshape(permute(AUZ2,[3,2,1]),[I3 R1*R2]);

ADUDVAW = (A3Z1V + A3UZ2)*AX.RW' + AX.A3UV*(A3Z1V' + A3UZ2')*W;

arg_of_skew = AX.invWTW*( ...
    -(Z3'*W + W'*Z3)*AX.invWTW_RWRW + Z3A3UVRW + Z3A3UVRW' + W'*ADUDVAW);

res{3} = Z3*AX.RWRW + W*(Z3A3UVRW + Z3A3UVRW') ...
    + W*(W'*ADUDVAW) - AX.A3UV_A3UV*Z3 -ADUDVAW ...
    + Z3*AX.skew_invWTW_RWRW + W*skew(arg_of_skew);
end
