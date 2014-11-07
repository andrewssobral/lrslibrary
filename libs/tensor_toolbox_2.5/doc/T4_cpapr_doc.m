%% Alternating Poisson Regression for fitting CP to sparse count data

%% Set up a sample problem
% We follow the general procedure outlined by E. C. Chi and T. G. Kolda, 
% On Tensors, Sparsity, and Nonnegative Factorizations, arXiv:1112.2414
% [math.NA], December 2011 (http://arxiv.org/abs/1112.2414).

% Pick the size and rank
sz = [100 80 60];
R = 5;

% Generate factor matrices with a few large entries in each column; this
% will be the basis of our soln.
A = cell(3,1);
for n = 1:length(sz)
    A{n} = rand(sz(n), R);
    for r = 1:R
        p = randperm(sz(n));
        nbig = round( (1/R)*sz(n) );
        A{n}(p(1:nbig),r) = 100 * A{n}(p(1:nbig),r);
    end
end
lambda = rand(R,1);
S = ktensor(lambda, A);
S = normalize(S,'sort',1);

% Create sparse test problem based on provided solution. 
nz = prod(sz) * .05;
info = create_problem('Soln', S, 'Sparse_Generation', nz);

% Extract data and solution
X = info.Data;
M_true = info.Soln;

%% Call CP-APR

% Compute a solution
M = cp_apr(X, R, 'printitn', 10);

% Score the solution
factor_match_score = score(M, M_true, 'greedy', true)
