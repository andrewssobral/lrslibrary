function [U,S,output] = lmlra_hooi(T,U0,options)
%LMLRA_HOOI LMLRA by higher-order orthogonal iteration.
%   [U,S,output] = lmlra_hooi(T,U0) computes the factor matrices U{1}, ...,
%   U{N} and core tensor S belonging to a low multilinear rank
%   approximation of the N-th order tensor T. The algorithm is initialized
%   with the factor matrices U0{n}. The structure output returns additional
%   information:
%
%      output.info       - The circumstances under which the procedure
%                          terminated:
%                             1: Objective function tolerance reached.
%                             2: Maximum number of iterations reached.
%      output.iterations - The number of iterations.
%      output.normS      - The Frobenius norm of the core tensor S in each
%                          iteration.
%      output.sangle     - The difference in subspace angle of U{1} between
%                          every two successive iterates.
%
%   lmlra_hooi(T,U0,options) may be used to set the following options:
%
%      options.MaxIter = 500 - The maximum number of iterations.
%      options.TolFun = 1e-6 - The tolerance for difference in subspace
%                              angle of U{1} between two successive
%                              iterates.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] P.M. Kroonenberg, J. de Leeuw, "Principal component analysis of
%       three-mode data by means of alternating least squares algorithms",
%       Psychometrika, Vol. 45, 1980, pp. 69-97.
%   [2] L. De Lathauwer, B. De Moor and J. Vandewalle, "On the best rank-1
%       and rank-(R1,R2,...,RN) approximation of higher-order tensors",
%       SIAM J. Matrix Anal. Appl., Vol. 21, 2000, pp. 1324-1342.
%   [3] P.M. Kroonenberg, "Applied Multiway Data Analysis", Wiley, 2008.

% Check the tensor T.
N = ndims(T);
if N < 3
    error('lmlra_hooi:T','ndims(T) should be >= 3.');
end

% Check the initial factor matrices U0.
U = U0(:).';
size_tens = size(T);
size_core = cellfun('size',U,2);
if any(cellfun('size',U,1) ~= size_tens)
    error('lmlra_hooi:U0','size(T,n) should equal size(U0{n},1).');
end

% Check the options structure.
if nargin < 3, options = struct; end
if ~isfield(options,'MaxIter'), options.MaxIter = 500; end
if ~isfield(options,'TolFun'), options.TolFun = 1e-6; end

% Higher-order orthogonal iteration.
output.info = false;
output.iterations = 0;
output.normS = [];
output.sangle = inf;
while ~output.info

    % Save current state.
    U1 = U{1};

    % Update factor matrices.
    for n = 1:N
        mode = [1:n-1 n+1:N];
        [T1,iperm] = tmprod(T,U(mode),mode,'H');
        T1 = reshape(permute(T1,iperm([n 1:n-1 n+1:N])),size_tens(n),[]);
        [Un,s,~] = svd(T1,'econ');
        size_core(n) = min(size_core(n),size(Un,2));
        U{n} = Un(:,1:size_core(n));
        s = diag(s); s = s(1:size_core(n));
    end

    % Update the output structure.
    output.iterations = output.iterations+1;
    output.normS(end+1) = norm(s);
    output.sangle(end+1) = subspace(U{1},U1);
    output.info = abs(output.sangle(end)) <= options.TolFun;
    if output.iterations >= options.MaxIter, output.info = 2; end

end

% Format output.
if nargout >= 2, S = mat2tens(U{N}'*T1,size_core,N); end
