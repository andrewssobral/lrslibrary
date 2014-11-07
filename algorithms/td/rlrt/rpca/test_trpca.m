function results = test_trpca( dataname, alg, rRatio, lambdaS, IsTC, verbose, mode, mu1fac, ks )
% dataname could be a string or actual data
% lambda1 = lambdaS * r * rRatio
% mode = RPCA mode (see rpca_for_tensor)
% mu1 = mu1fac * std(T(:))

if ~isstruct( dataname )
    load( ['..\..\data\', dataname] );
else
    data = dataname;
end

algs = get_algs();
params.X0 = tenzeros( size(data.T) );
N = ndims(data.T);
params.V0 = cell( 1, N );
for i = 1:N
    params.V0{i} = tenzeros( size(data.T) );
end
params.E0 = tenzeros( size(data.T) );
params.mu0 = 1/(N+1);
T = double(data.T);
if ~exist( 'mu1fac', 'var' ) || isempty(mu1fac); mu1fac = 10; end
params.mu1fac = mu1fac;
params.mu1 = mu1fac*std(data.b);  %mu1fac*std(T(:));
params.mu2 = params.mu1;
params.mu_min = 1e-4;
params.mu_max = 1e2;
params.max_iter = 1000;
params.opt_tol = 1e-3;  %1e-5 for syn data analysis plots, 1e-3 for others
params.eta = 1/(N+1);
params.IsTC = exist( 'IsTC', 'var' ) && ~isempty(IsTC) && IsTC;

r = 1/sqrt( max(size(data.T)) ); %0.015; %0.09:-0.005:0.02;
if ~exist( 'lambdaS', 'var' ) || isempty(lambdaS)
    lambdaS = 1;
end
params.lambdaS = lambdaS;     %1e1;
if ~exist( 'rRatio', 'var' ) || isempty(rRatio)
    rRatio = 1/4;
end
params.lambda = params.lambdaS*r*rRatio;       %params.lambdaS*r/4;
% params.lambda = 0.06;
params.rRatio = rRatio;
params.verbose = exist( 'verbose', 'var' ) && ~isempty(verbose) && verbose;
params.use_cont = true;
if ~exist( 'mode', 'var' ) || isempty(mode)
    mode = N;
end

%%%%%%%%%% for PROPACK %%%%%%%%%%%%
% declare global var 'sv'
global sv;
global tmode;
global use_propack;
global curr_mu;
sv =  ceil(min(size(data.T)) * 0.1) * ones( 1, N );
use_propack = true;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch alg
    case 1
        results = tensor_rpca_adal2( data, params );
    case 2
        results = tensor_rpca_tc_adal( data, params );
    case 3
        params.mode = mode;
        results = rpca_for_tensor( data, params );
   
end

results.rRatio = rRatio;
results.alg = algs{alg};
results.dataname = dataname;
results.X0 = data.X;
if isfield( data, 'noise' )
    results.noise = data.noise;
end
[ rel_err, snr ] = get_SNR( results.X, results.X0 );
results.rel_err1 = rel_err;
results.snr1 = snr;


% save hall_results results;

% % run matrix version of RPCA on vectorized version of data
% if exist( 'runRPCA', 'var' ) && ~isempty(runRPCA) && runRPCA
%     results = rpca_for_tensor( data, params, results );
%     [ rel_err, snr ] = get_SNR( results.A, data.X );
%     results.rel_err2 = rel_err;
%     results.snr2 = snr;
% end



end
