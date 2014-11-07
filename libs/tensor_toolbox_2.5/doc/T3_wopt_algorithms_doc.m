%% Weighted optimization for CP tensor decomposition with incomplete data
% We explain how to use |cp_wopt| with the POBLANO toolbox. 

%% Create an example problem with missing data. 
% Here we have 25% missing data and 10% noise.   
R = 2;
info = create_problem('Size', [15 10 5], 'Num_Factors', R, ...
    'M', 0.25, 'Noise', 0.10);
X = info.Data;
P = info.Pattern;
M_true= info.Soln;

%% Create initial guess using 'nvecs'
M_init = create_guess('Data', X, 'Num_Factors', R, ...
    'Factor_Generator', 'nvecs');


%% Set up the optimization parameters
% It's genearlly a good idea to consider the parameters of the optimization
% method. The default options may be either too stringent or not stringent
% enough. The most important options to consider are detailed here. 

% Get the defaults
ncg_opts = ncg('defaults');
% Tighten the stop tolerance (norm of gradient). This is often too large.
ncg_opts.StopTol = 1.0e-6;
% Tighten relative change in function value tolearnce. This is often too large.
ncg_opts.RelFuncTol = 1.0e-20;
% Increase the number of iterations. 
ncg_opts.MaxIters = 10^4;
% Only display every 10th iteration
ncg_opts.DisplayIters = 10;
% Display the final set of options
ncg_opts

%% Call the |cp_wopt| method
% Here is an example call to the cp_opt method. By default, each iteration
% prints the least squares fit function value (being minimized) and the
% norm of the gradient. The meaning of any line search warnings
% can be checked via <matlab:doc('cvsrch') doc cvsrch>.
[M,~,output] = cp_wopt(X, P, R, 'init', M_init, ...
    'alg', 'ncg', 'alg_options', ncg_opts);

%% Check the output
% It's important to check the output of the optimization method. In
% particular, it's worthwhile to check the exit flag. 
% A zero (0) indicates successful termination with the gradient smaller
% than the specified StopTol, and a three (3) indicates a successful
% termination where the change in function value is less than RelFuncTol.
% The meaning of any other flags can be checked via 
% <matlab:doc('poblano_params') doc poblano_params>. 
exitflag = output.ExitFlag


%% Evaluate the output
% We can "score" the similarity of the model computed by CP and compare
% that with the truth. The |score| function on ktensor's gives a score in
% [0,1]  with 1 indicating a perfect match. Because we have noise, we do
% not expect the fit to be perfect. See <matlab:doc('ktensor/score') doc
% score> for more details.
scr = score(M,M_true)

%% Create a SPARSE example problem with missing data. 
% Here we have 95% missing data and 10% noise.   
R = 2;
info = create_problem('Size', [150 100 50], 'Num_Factors', R, ...
    'M', 0.95, 'Sparse_M', true, 'Noise', 0.10);
X = info.Data;
P = info.Pattern;
M_true= info.Soln;

%% Create initial guess using 'nvecs'
M_init = create_guess('Data', X, 'Num_Factors', R, ...
    'Factor_Generator', 'nvecs');


%% Set up the optimization parameters
% It's genearlly a good idea to consider the parameters of the optimization
% method. The default options may be either too stringent or not stringent
% enough. The most important options to consider are detailed here. 

% Get the defaults
ncg_opts = ncg('defaults');
% Tighten the stop tolerance (norm of gradient). This is often too large.
ncg_opts.StopTol = 1.0e-6;
% Tighten relative change in function value tolearnce. This is often too large.
ncg_opts.RelFuncTol = 1.0e-20;
% Increase the number of iterations. 
ncg_opts.MaxIters = 10^4;
% Only display every 10th iteration
ncg_opts.DisplayIters = 10;
% Display the final set of options
ncg_opts

%% Call the |cp_wopt| method
% Here is an example call to the cp_opt method. By default, each iteration
% prints the least squares fit function value (being minimized) and the
% norm of the gradient. The meaning of any line search warnings
% can be checked via <matlab:doc('cvsrch') doc cvsrch>.
[M,~,output] = cp_wopt(X, P, R, 'init', M_init, ...
    'alg', 'ncg', 'alg_options', ncg_opts);

%% Check the output
% It's important to check the output of the optimization method. In
% particular, it's worthwhile to check the exit flag. 
% A zero (0) indicates successful termination with the gradient smaller
% than the specified StopTol, and a three (3) indicates a successful
% termination where the change in function value is less than RelFuncTol.
% The meaning of any other flags can be checked via 
% <matlab:doc('poblano_params') doc poblano_params>. 
exitflag = output.ExitFlag


%% Evaluate the output
% We can "score" the similarity of the model computed by CP and compare
% that with the truth. The |score| function on ktensor's gives a score in
% [0,1]  with 1 indicating a perfect match. Because we have noise, we do
% not expect the fit to be perfect. See <matlab:doc('ktensor/score') doc
% score> for more details.
scr = score(M,M_true)

