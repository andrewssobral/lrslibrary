function [U,output] = cpd3_sd(T,U0,options)
%CPD3_SD CPD by simultaneous diagonalization.
%   [U,output] = cpd3_sd(T,U0) computes the factor matrices U{1}, U{2} and
%   U{3} belonging to a canonical polyadic decomposition of the third-order
%   tensor T. The algorithm is initialized with the factor matrices U0{n}.
%   The number of rank-one terms R should satisfy 1 <= R <= max(size(T)).
%   The structure output returns additional information:
%
%      output.SDName   - The name of the selected simultaneous
%                        diagonalization algorithm.
%      output.SDOutput - Information returned by the selected simultaneous
%                        diagonalization algorithm.
%
%   cpd3_sd(T,U0,options) may be used to set the following options:
%
%      options.SD =            - A function handle to the desired
%      [@cpd_als|@cpd_minf|...   simultaneous diagonalization algorithm.
%       @cpd_nls|{@cpd3_sgsd}]
%      options.SDOptions       - An options structure passed to the
%                                selected SD algorithm. See also
%                                help [options.SD].

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] L. De Lathauwer, "A Link between the Canonical Decomposition in 
%       Multilinear Algebra and Simultaneous Matrix Diagonalization,"
%       SIAM J. Matrix Anal. Appl., Vol. 28, 2006, pp. 642-666.

% Incomplete/sparse tensors not supported yet.
T = fmt(T,true);
if isstruct(T)
    error('This method does not yet support incomplete/sparse tensors.');
end

% Check the tensor T.
N = ndims(T);
if N ~= 3, error('cpd3_sd:T','ndims(T) should be 3.'); end

% Check the initial factor matrices U0.
U0 = U0(:).';
R = size(U0{1},2);
size_tens = size(T);
if any(cellfun('size',U0,2) ~= R)
    error('cpd3_sd:U0','size(U0{n},2) should be the same for all n.');
end
if R > max(size_tens)
    error('cpd3_sd:U0','size(U0{n},2) should be <= max(size(T)).');
end
if any(cellfun('size',U0,1) ~= size_tens)
    error('cpd3_sd:U0','size(T,n) should equal size(U0{n},1).');
end

% Check the options structure.
xsfunc = @(f)isa(f,'function_handle')&&exist(func2str(f),'file');
if nargin < 3, options = struct; end
if ~isfield(options,'SD')
    options.SD = @cpd3_sgsd;
elseif ~xsfunc(options.SD)
    error('cpd3_sd:SD','Algorithm does not exist.');
end

% Prepermute the tensor so the longest mode is the first one, and for
% efficiency, the last one is the shortest one.
[size_tens,idx] = sort(size_tens,'descend');
if any(idx ~= 1:3)
    T = permute(T,idx);
    U0 = U0(idx);
    iperm(idx) = 1:N;
end
    
% Best rank-R approximation of the mode-1 unfolding of T.
% T(1) = U{1}*kr(U([3 2]))^T = F*W^-H * W^H*E^H
% U{1} = F*W^-H; kr(U([3 2])) = conj(E*W)
T1 = reshape(T,size_tens(1),[]);
[~,~,E] = svd(T1,'econ');
E = E(:,1:R);

% Compute the Gramian P'*P for the rank-one detecting device P. The Gramian
% is computed as a fourth-order tensor in which element with indices
% (r,s,t,u) is equal to <Er,Et>*<Es,Eu> + <Er,Eu>*<Es,Et> - <Eu'*Es,Er'*Et>
% - <Et'*Es,Er'*Eu>, where En denotes the matricized n-th column of E.
% First compute all products P34(r,s,t,u) = conj(<Er'*Es,Et'*Eu>).
EE = reshape(E,size_tens(2),[]);
P34 = zeros(size_tens(3)^2,R^2);
for r = 1:R
    X = EE(:,(r-1)*size_tens(3)+(1:size_tens(3)))'*EE;
    P34(:,(r-1)*R+(1:R)) = reshape(X,size_tens(3)^2,R);
end
P34 = reshape(P34'*P34,[R R R R]);
% Compute all four terms P1 to P4.
rr = false(1,R*(R+1)/2);
rr([1 1+cumsum(R:-1:2)]) = true;
P12 = ones(R*(R+1)/2,1);
P12(rr) = 2;
P12 = diag(P12);
P34 = -reshape(permute(P34,[1 4 2 3]),[R^2 R^2]) ...
      -reshape(permute(P34,[1 4 3 2]),[R^2 R^2]);
% Compose the matrix P'*P, without redundant columns/rows.
lt = tril(true(R,R));
PP = P12+P34(lt(:),lt(:));

% Compute the kernel of the adjusted Gramian and build the tensor M,
% consisting of symmetric matrices Mk = W Lk W^T.
[~,~,K] = svd(PP,'econ');
K = K(:,end:-1:end-R+1);
K(~rr,:) = 0.5*K(~rr,:);
M = zeros(R,R,R);
Mlt = zeros(R,R);
for r = 1:R
    Mlt(lt) = K(:,r);
    M(:,:,r) = Mlt+tril(Mlt,-1).';
end

% Simultaneously diagonalize the frontal slices of M.
W0 = E\conj(kr(U0(end:-1:2)));
C0 = reshape(M,[R*R R]).'*conj(kr(W0,W0))/conj((W0'*W0).^2);
if ~isfield(options,'SDOptions'), options.SDOptions = struct; end
[W,output.SDOutput] = options.SD(M,{W0,W0,C0},options.SDOptions);
output.SDName = func2str(options.SD);

% Reconstruct the factor matrices.
U = cell(1,N);
krCB = conj(E*W{1});
for r = 1:R
    [u,s,v] = svd(reshape(krCB(:,r),size_tens([2 3])),'econ');
    U{2} = [U{2} s(1)*u(:,1)];
    U{3} = [U{3} conj(v(:,1))];
end
U{1} = T1*conj(kr(U([3 2])))/conj((U{3}'*U{3}).*(U{2}'*U{2}));

% Inverse permute the factor matrices.
if exist('iperm','var')
    U = U(iperm);
end
