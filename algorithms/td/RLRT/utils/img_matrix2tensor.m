function T = img_matrix2tensor( M, n1, n2, maxPix )

N = size( M, 2 );
T = zeros( n1, n2, N );

for i = 1:N
    T(:,:,i ) = reshape( M(:,i), [ n1, n2 ] ) / maxPix;
end

T = tensor( T );

end