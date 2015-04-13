function [U,output] = cpd_gevd(T,R,options)
%CPD_GEVD CPD by a generalized eigenvalue decomposition.
%   [U,output] = cpd_gevd(T,R) computes the factor matrices U{1}, ..., U{N}
%   belonging to a canonical polyadic decomposition of the N-th order
%   tensor T. The number of rank-one terms R should satisfy
%   sum(R <= size(T)) >= 2. The decomposition is exact if T is rank-R and
%   two factor matrices have full column rank. If this is not the case, the
%   method may fall back on a different initialization, see options.Backup.
%
%   cpd_gevd(T,R,options) may be used to set the following options:
%
%      options.Backup =      - The method assumes the first two factor
%      @cpd_rnd                matrices have full column rank. If the
%                              method detects this is not the case, it will
%                              fall back to the initialization
%                              options.Backup(T,R,options).
%      options.IsReal =      - If true, the factor matrices U{n} are forced
%      [{'auto'}|false|true]   to be real. If false, complex factor
%                              matrices are allowed. On 'auto',
%                              options.IsReal is equal to isreal(T).

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)
%
%   References:
%   [1] S.E. Leurgans, R.T. Ross, R.B. Abel, "A decomposition for three-way
%       arrays," SIAM J. Matrix Anal. Appl., Vol. 14, 1993, pp. 1064–1083.
%   [2] E. Sanchez, B.R. Kowalski, "Tensorial resolution: A direct
%       trilinear decomposition," J. Chemometrics, Vol. 4, 1990, pp. 29–45.
%   [3] R. Sands, F. Young, "Component models for three-way data: An
%       alternating least squares algorithm with optimal scaling features,"
%       Psychometrika, Vol. 45, 1980, pp. 39–67.

% Incomplete/sparse tensors not supported yet.
T = fmt(T,true);
if isstruct(T)
    error('This method does not yet support incomplete/sparse tensors.');
end

% Check the tensor T.
N = ndims(T);
if N < 3, error('cpd_gevd:T','ndims(T) should be >= 3.'); end

% Check the number of rank one terms R.
size_tens = size(T);
if sum(R > size_tens) > 1
    error('cpd_gevd:R','sum(size(U0{n},2) <= size(T)) should be >= 2.');
end

% Check the options structure.
if nargin < 3, options = struct; end
if ~isfield(options,'Backup'), options.Backup = @cpd_rnd; end
if ~isfield(options,'IsReal'), options.IsReal = 'auto'; end
if ischar(options.IsReal) && strcmpi(options.IsReal,'auto')
    options.IsReal = isreal(T);
end

% Compute a direct multilinear decomposition.
if R > 1
    
    % Prepermute the (core) tensor so the largest two modes are the first
    % two, which hopefully maximizes the probability that the first two
    % factor matrices have full column rank.
    [V,S,sv] = mlsvd(T);
    tol = max(size(T),numel(T)./size(T))*eps(class(T));
    size_core = arrayfun(@(n)sum(sv{n} > tol(n)*max(sv{n})),1:ndims(T));
    [size_core,idx] = sort(size_core,'descend');
    T = permute(T,idx);
    size_tens = size_tens(idx);
    S = permute(S,idx);
    V = V(idx);
    iperm(idx) = 1:N;
    
    % Test to see if the method can be applied.
    size_core(1:2) = min(size_core(1:2),R);
    size_core(3) = min(size_core(3),2);
    size_core(4:end) = 1;
    idx = arrayfun(@(i)1:i,size_core,'UniformOutput',false);
    S = S(idx{:});
    for n = 1:N, Vn = V{n}; V{n} = Vn(:,1:size_core(n)); end
    if any(R > size_core(1:2)) || size_core(3) ~= 2
        % Assumption of two full column rank matrices can not be met,
        % fall back to a random initialization.
        [U,output] = options.Backup(T,R,options);
        output.Name = func2str(options.Backup);
        return;
    end
    
    % Retrieve U{1} using a GEVD.
    [Apt,Dpt] = eig(S(:,:,1).',S(:,:,2).');
    if options.IsReal && ~isreal(Apt)
        c = find(imag(diag(Dpt)) ~= 0);
        for i = 1:2:length(c)
            Apt(:,c(i)) = real(Apt(:,c(i)));
            Apt(:,c(i+1)) = imag(Apt(:,c(i+1)));
        end
    end
    
    % Retrieve remaining factor matrices using the MLSVD.
    U = cell(1,N);
    U{1} = V{1}/(Apt.');
    T1 = reshape(T,size_tens(1),[]);
    X = T1.'*conj(V{1})*Apt;
    for r = 1:R
        [u,s] = mlsvd(reshape(X(:,r),size_tens(2:end)),ones(1,N-1));
        u{1} = u{1}*s(1);
        for n = 2:N, U{n} = [U{n} u{n-1}]; end
    end
    U{1} = T1/(kr(U(end:-1:2)).');
    
else
    
    % Compute an approximate best rank-one approximation.
    [U,S] = mlsvd(T,ones(1,N));
    U{1} = U{1}*S;
    
end

% Inverse permute the factor matrices.
if exist('iperm','var')
    U = U(iperm);
end
output.Name = 'cpd_gevd';
