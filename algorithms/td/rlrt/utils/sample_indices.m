function [ b, Omega, linInd ] = sample_indices( X, frac )

dims = size( X );
Ndims = length(dims);
N = 1;
for i = 1:Ndims
    N = N * dims(i);
end

Nfrac = ceil( N * frac );
Irand = randperm( N )';
linInd = Irand( 1:Nfrac );

if Ndims == 3
    [I,J,K] = ind2sub( dims, linInd );
    Omega = [ I J K ];
else if Ndims == 2
        [I,J] = ind2sub( dims, linInd );
        Omega = [ I J ];
    end
end
b = X( Omega );

end