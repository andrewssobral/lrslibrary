function data = add_batch_noise_obs( dataname, obsPct, noisePct, mag, Nrep, data )
% obsPct = 5%, 10%, ..., 95%, 100%
% noisePct = 1%, 5%, 10%, ..., 45%
% mag = 1, 5
%
% Data structure:
% Observations are stored as vectors of linear indices.
% Noise values are stored as (value vector, index vector) pairs.
% Indices are linear indices.
% Each observation and noise % value has 5 random replicates.
% obs{i} = Nt x Nrep matrix of indices
% noise{i} = vector of length Nrep
% noise{i}(k).val(:,j) = value vector of kth replicate with j-th magnitude
% noise{i}(k).ind = index vector of kth replicate

% generate original tensor if not given in input
if ~exist( 'data', 'var' ) || isempty(data)
    data = gen_lowrank_tensor( [50 50 20], [5 5 5] );
end

X = double(data.X);
N = length(X(:));
Nrep = 5;

% set magnitudes of noise
mag = [ 1, 5 ];
Nmag = length(mag);

% sample observations
obsPct = 0.05 : 0.05 : 1;
Nobs = length(obsPct);
obs = cell( 1, Nobs );
for i = 1:Nobs
    Nt = ceil( N*obsPct(i) );
    obs{i} = zeros( Nt, Nrep );
    for k = 1:Nrep
        Irand = randperm(N)';
        obs{i}(:,k) = Irand(1:Nt);
    end
end

% generate noise
noisePct = [ 0.01, 0.05:0.05:0.45 ];
Nns = length(noisePct);
noise = cell( 1, Nns );
for i = 1:Nns
    Nt = ceil( N*noisePct(i) );
    noise{i} = [];
    for k = 1:Nrep
        noise{i}(k).val = zeros( Nt, Nmag );
        for j = 1:Nmag
            noise{i}(k).val(:,j) = random( 'unif', -mag(j), mag(j), Nt, 1 );
        end
        Irand = randperm(N)';
        noise{i}(k).ind = Irand(1:Nt);
    end
end

data.mag = mag;
data.obsPct = obsPct;
data.obs = obs;
data.noisePct = noisePct;
data.noise = noise;
data.Nrep = Nrep;

save( ['..\..\data\', dataname], 'data' );
end