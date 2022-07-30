%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Rank Sparsity Tensor Decomposition 
% Reference: "Optimum Subspace Learning and Error Correction for Tensors",ECCV 2010
%
% min(L, S, M1,...Mn, N1,...Nn): 
%    (alpha1||L_(1)-M1||^2 + alpha2||L_(2)-M2||^2 + alpha3||X_(3)-M3||^2 + ...)/(2n) + 
%    (beta1||S_(1)-N1||^2 + beta2||S_(2)-N2||^2 + beta3||S_(3)-N3||^2 + ...)/(2n) +
%    (gamma1||X_(1)-N1-M1||^2 + gamma2||X_(2)-N2-M2||^2 + gamma3||X_(3)-N3-M3||^2 + ...)/(2n) +
%    (lambda1||M1||_tr + lambda2||M2||_tr + lambda3||M3||_tr +)/n ... +
%    (eta||N1||_1 + eta||N2||_1 + eta||N3||_1 + ...)/n
%
% Version:   0.9
% Revised:   05/04/2009
% Reference: Optimum Subspace Learning and Error Correction for Tensors, ECCV 2010
% Author:    Yin Li @ SJTU
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L, S, Ud, rank, sparsity, errorList, iter] = RSTD(varargin)
%%%
% parameter checking
X = varargin{1};

alpha = varargin{2};    % ||M_i-L_{(i)}||_F
beta = varargin{3};     % ||N_i-S_{(i)}||_F
gamma = varargin{4};    % ||M_i+N_i-X_{(i)}||_F
 
lambda = varargin{5};   % ||M||_tr
eta = varargin{6};      % ||N||_1

maxIter = varargin{7};

sv = varargin{8};

tol = 1e-6;

% we do not check the parameter in current version

%%%
%initialization
dim = size(X);
M = cell(ndims(X), 1); %sets of matrix
N = cell(ndims(X), 1); %sets of matrix
Ud = cell(ndims(X),1);
n_dim = zeros([ndims(X) 1]);
svp = sv;

for i = 1:ndims(X)
    M{i} = unfolding(X,i);
    N{i} = zeros(size(M{i}));
    n_dim(i) = min(size(M{i}));
end
L = X;
S = zeros(dim);

rank = zeros([ndims(X) maxIter]);
sparsity = zeros([ndims(X) maxIter]);
errorList = zeros([2 maxIter]);

Lsum = zeros(dim);
Ssum = zeros(dim);

ACsum = alpha + gamma;
BCsum = beta + gamma;

t = 1;

%%%
for iter = 1:maxIter
    %iter
     
    %clear
    OldL = L;
    OldS = S;
    Lsum = Lsum .* 0;
    Ssum = Ssum .* 0;
    
    for i = 1:ndims(X)
        
        Mpro = (alpha(i) .* unfolding(L,i) + gamma(i) .* (unfolding(X,i)-N{i}))./ACsum(i);
        Npro = ( beta(i) .* unfolding(S,i) + gamma(i) .* (unfolding(X,i)-M{i}))./BCsum(i);
                 
        % subproblem for N{i} using l1 norm 
        N{i} = epp1(Npro, eta(i)/BCsum(i));
        sparsity(i,iter) = (length(find(abs(N{i})>tol)))./numel(Npro);
                           
        % subproblem for M{i} using trace norm
        [U, sigma, V] = smartSVD(Mpro, sv(i));
        sigma = epp1(diag(sigma), lambda(i)/ACsum(i)); %sigma = max(diag(sigma) - lambda(i)/ACsum(i), 0);
        svp(i) = sum(sigma > tol);
        M{i} = U(:, 1:svp(i)) * diag(sigma(1:svp(i))) * V(:, 1:svp(i))';
        Ud{i} = U(:, 1:svp(i));
        rank(i,iter) = svp(i);
        
        % predict the dimension of the subspace
        if svp(i) < sv(i)
            sv(i) = min(svp(i) + 1, n_dim(i));
        else
            sv(i) = min(sv(i) + round(0.05*n_dim(i)), n_dim(i));   
        end
        
        Ssum = Ssum +  beta(i) .* folding(N{i}, i, dim);
        Lsum = Lsum + alpha(i) .* folding(M{i}, i, dim);
               
    end
    
    S = Ssum ./ (sum(beta) + 1e-10);
    L = Lsum ./ (sum(alpha) + 1e-10);
    
    if iter>1
        Serror = sum(sum(unfolding(abs((OldS-S)),1)))./numel(S);
        Lerror = sum(sum(unfolding(abs((OldL-L)),1)))./numel(L);
        errorList(1,iter) = Lerror;
        errorList(2,iter) = Serror;
        
        if Serror <= tol && Lerror<=tol
            disp(['Converge at ', num2str(iter), ' iterations with the tolerance of ', num2str(tol)]);
            rank = rank(:,1:iter);
            sparsity = sparsity(:,1:iter);
            errorList = errorList(:,1:iter);
            return;
        end  
        
    end
    
    % implementation of proximal gradient, it is much faster (upto 40%) than the original BCD based method
    Lsum = Lsum .* 0;
    Ssum = Ssum .* 0;
    t_new = 0.5.*(1+sqrt(1+4*t^2));
    for i=1:ndims(X)         
        N{i} = unfolding(S,i) + (1./t_new).*(t.*N{i}+(1-t).*unfolding(OldS,i)-unfolding(S,i));
        M{i} = unfolding(L,i) + (1./t_new).*(t.*M{i}+(1-t).*unfolding(OldL,i)-unfolding(L,i));
        Ssum = Ssum +  beta(i) .* folding(N{i}, i, dim);
        Lsum = Lsum + alpha(i) .* folding(M{i}, i, dim);
    end
    S = Ssum ./ (sum(beta) + 1e-10);
    L = Lsum ./ (sum(alpha) + 1e-10);
    t = t_new; 
    % end proximal gradient
    
end

disp(['Reaches ', num2str(iter), 'iterations with the tolerance of', num2str(tol)]);
rank = rank(:,1:iter);
sparsity = sparsity(:,1:iter);
errorList = errorList(:,1:iter);
