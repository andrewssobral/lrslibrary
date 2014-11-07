function data = add_sparse_gross_noise( data_name, spa, mag, frac, doSave )
% add sparse impulsive noise
% sample frac portion of data

if isstruct( data_name )
    data = data_name;
    doSave = false;
else
    load( ['..\..\data\',data_name] );
end
N = length( size(data.X) );
x = tenmat( data.X, [1:N] );
m = size( x, 1 );
I = randperm( m )';
Nnoise = ceil(m*spa);
noiseI = I( 1:Nnoise );
x( noiseI ) = x( noiseI ) + random( 'unif', -mag, mag, Nnoise, 1 );
data.T = tensor( x );
data.noiseI = noiseI;
data.noise = mag;

Nobs = ceil(m*frac);
linInd = I( 1:Nobs );
data.b = data.T(linInd);
data.linInd = linInd;
data.frac = frac;

if ~exist( 'doSave', 'var' ) || doSave
    save( ['..\..\data\',data_name,'-n'], 'data' );
end
end