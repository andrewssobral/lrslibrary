function results = test_tc( dataname )

load( ['..\..\data\', dataname] );
addpath( '..\utils' );

% data.W = tenzeros( size(data.X) );
% data.W(data.Omega) = data.b;
data.W = tenzeros( size(data.X) );
N = length( size(data.W) );
data.Lamb = cell(1,N);
for i = 1:N
    data.Lamb{i} = tenzeros( size(data.X) );
end
data.Lambda = tenzeros( size(data.X) );

params.mu0 = 1e1;
params.sigma = 1e-3;
params.lambda = 1;
params.max_iter = 100;
params.opt_tol = 1e-3;
params.R = 50;

while params.lambda < 100
    results = TC_ADAL_Nuc( data, params );
%     results = TC_MSA( data, params );
    data.W = results.X;
    params.lambda = params.lambda * 5;
end

end