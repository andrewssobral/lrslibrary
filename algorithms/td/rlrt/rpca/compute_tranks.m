function [ rankSum, normRank ] = compute_tranks( Rs, rs )

N = length(rs);
temp = zeros( 1, N );
for i = 1:N
    temp(i) = min( rs(i), prod(rs)/rs(i) );
end
rankSum = sum( temp );
normRank = halfNorm( 1./Rs ) * halfNorm( rs );

end

function val = halfNorm( x )

val = mean( sqrt(x) )^2;

end