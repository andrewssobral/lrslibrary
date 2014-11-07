%% Creating test problems and initial guesses
% We demonstrate how to use Tensor Toolbox |create_problem| and
% |create_guess| functions to create test problems for fitting algorithms. 

%% Creating a CP test problem
% The |create_problem| function allows a user to generate a test problem
% with a known solution having a pre-specified solution. The
% |create_problem| function generates both the solution (as a |ktensor| for
% CP) and the test data (as a |tensor|). We later show that a
% pre-specificed solution can be used as well.

% Create a problem
info = create_problem('Size', [5 4 3], 'Num_Factors', 3, 'Noise', 0.10);

%%

% Display the solution created by create_problem
soln = info.Soln

%%

% Display the data created by create_problem
data = info.Data

%%

% The difference between true solution and measured data should match the
% specified 10% noise.
diff = norm(full(info.Soln) - info.Data)/norm(full(info.Soln))

%% Creating a Tucker test problem
% The |create_problem| function can also be used to create Tucker problems
% by specifying the |'Type'| as |'Tucker'|. In this case, the
% |create_problem| function generates both the solution (as a |ttensor| for
% Tucker) and the test data (as a |tensor|). 

% Create a problem
info = create_problem('Type', 'Tucker', 'Size', [5 4 3], 'Num_Factors', [3 3 2]);

%%

% Display the Tucker-type solution created by create_problem
soln = info.Soln

%%

% Difference between true solution and measured data (default noise is 10%)
diff = norm(full(info.Soln) - info.Data)/norm(full(info.Soln))

%% Recreating the same test problem
% We can recreate exactly the same test problem when we use the same random
% seed and other parameters.

% Set-up, including specifying random state
sz = [5 4 3]; %<- Size
nf = 2; %<- Number of components
state = RandStream.getDefaultStream.State; %<- Random state

%%

% Generate first test problem
info1 = create_problem('Size', sz, 'Num_Factors', nf, 'State', state);

%%

% Generate second identical test problem
info2 = create_problem('Size', sz, 'Num_Factors', nf, 'State', state);

%%

% Check that the solutions are identical
tf = isequal(info1.Soln, info2.Soln)

%%

% Check that the data are identical
diff = norm(info1.Data - info2.Data)

%% Checking default parameters and recreating the same test problem
% The |create_problem| function returns the parameters that were used to
% generate it. These can be used to see the defaults. Additionally, if
% these are saved, they can be used to recreate the same test problems for
% future experiments.

% Generate test problem and use second output argument for parameters.
[info1,params] = create_problem('Size', [5 4 3], 'Num_Factors', 2);

%%

% Here are the parameters
params 

%%

% Recreate an identical test problem
info2 = create_problem(params);

%%

% Check that the solutions are identical
tf = isequal(info1.Soln, info2.Soln)

%%

% Check that the data are identical
diff = norm(info1.Data - info2.Data)

%% Options for creating factor matrices, core tensors, and lambdas
% Any function with two arguments specifying the size can be used to
% generate the factor matrices. This is specified by the
% |'Factor_Generator'| option to |create_problem|.
%
% Pre-defined options for |'Factor_Generator'| for creating factor matrices
% (for CP or Tucker) include:  
%
% * |'rand'| - Uniform on [0,1] 
% * |'randn'| - Gaussian with mean 0 and std 1
% * |'orthogonal'| - Generates a random orthogonal matrix. This option only
% works when the number of factors is less than or equal to the smallest
% dimension.
% * |'stochastic'| - Generates nonnegative factor matrices so that each
% column sums to one. 
%
% Pre-defined options for |'Lambda_Generator'| for creating lambda vector
% (for CP) include: 
%
% * |'rand'| - Uniform on [0,1] 
% * |'randn'| - Gaussian with mean 0 and std 1
% * |'orthogonal'| - Creates a random vector with norm one.
% * |'stochastic'| - Creates a random nonnegative vector whose entries sum
% to one. 
%
% Pre-defined options for |'Core_Generator'| for creating core tensors (for
% Tucker) include: 
%
% * |'rand'| - Uniform on [0,1] 
% * |'randn'| - Gaussian with mean 0 and std 1

% Here is ane example of a custom factor generator 
factor_generator = @(m,n) 100*rand(m,n);
info = create_problem('Size', [5 4 3], 'Num_Factors', 2, ...
    'Factor_Generator', factor_generator, 'Lambda_Generator', @ones);
first_factor_matrix = info.Soln.U{1}

%%

% Here is an example of a custom core generator for Tucker:
info = create_problem('Type', 'Tucker', 'Size', [5 4 3], ...
    'Num_Factors', [2 2 2], 'Core_Generator', @tenones);
core = info.Soln.core

%%

% Here's another example for CP, this time using a function to create
% factor matrices such that the inner products of the columns are
% prespecified.
info = create_problem('Size', [5 4 3], 'Num_Factors', 3, ...
    'Factor_Generator', @(m,n) tt_ccong(m,n,.9));
U = info.Soln.U{1};
congruences = U'*U

%% Generating data from an existing solution
% It's possible to skip the solution generation altogether and instead just
% generate appropriate test data.

% Manually generate a test problem (or it comes from some
% previous call to |create_problem|.
soln = ktensor({rand(50,3), rand(40,3), rand(30,3)});

% Use that soln to create new test problem.
info = create_problem('Soln', soln);

% Check whether solutions is equivalent to the input
iseq = isequal(soln,info.Soln)

%% Creating dense missing data problems
% It's possible to create problems that have a percentage of missing data.
% The problem generator randomly creates the pattern of missing data.

% Specify 25% missing data as follows:
[info,params] = create_problem('Size', [5 4 3], 'M', 0.25);

%% 

% Here is the pattern of known data (1 = known, 0 = unknown)
info.Pattern

%%

% Here is the data (incl. noise) with missing entries zeroed out
info.Data 

%% Creating sparse missing data problems. 
% If |Sparse_M| is set to true, then the data returned
% is sparse. Moreover, the dense versions are never explicitly created.
% This option only works when M >= 0.8.

% Specify 80% missing data and sparse
info = create_problem('Size', [5 4 3], 'M', 0.80, 'Sparse_M', true);

%% 

% Here is the pattern of known data
info.Pattern

%%

% Here is the data (incl. noise) with missing entries zeroed out
info.Data 

%% Create missing data problems with a pre-specified pattern
% It's also possible to provide a specific pattern (dense or sparse) to be
% used to specify where data should be missing.

% Create pattern
P = tenrand([5 4 3]) > 0.5;
% Create test problem with that pattern
info = create_problem('Size', size(P), 'M', P);
% Show the data
info.Data

%% Creating sparse problems (CP only)
% If we assume each model parameter is the input to a Poisson process, then
% we can generate a sparse test problems. This requires that all the factor
% matrices and lambda be nonnegative. The default factor generator
% ('randn') won't work since it produces both positive and negative values.

% Generate factor matrices with a few large entries in each column; this
% will be the basis of our soln.
sz = [20 15 10];
nf = 4;
A = cell(3,1);
for n = 1:length(sz)
    A{n} = rand(sz(n), nf);
    for r = 1:nf
        p = randperm(sz(n));
        idx = p(1:round(.2*sz(n)));
        A{n}(idx,r) = 10 * A{n}(idx,r);
    end
end
S = ktensor(A);
S = normalize(S,'sort',1);
%%

% Create sparse test problem based on provided solution. The
% 'Sparse_Generation' says how many insertions to make based on the
% provided solution S. The lambda vector of the solution is automatically
% rescaled to match the number of insertions.
info = create_problem('Soln', S, 'Sparse_Generation', 500);
num_nonzeros = nnz(info.Data)
total_insertions = sum(info.Data.vals)
orig_lambda_vs_rescaled = S.lambda ./ info.Soln.lambda

%% Generating an initial guess
% The |create_guess| function creates a random initial guess as a cell
% array of matrices. Its behavior is very similar to |create_problem|. A
% nice option is that you can generate an initial guess that is a
% pertubation of the solution.

info = create_problem;

% Create an initial guess to go with the problem that is just a 5%
% pertubation of the correct solution.
U = create_guess('Soln', info.Soln, 'Factor_Generator', 'pertubation', ...
    'Pertubation', 0.05);

