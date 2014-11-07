function data = prep_tensor( impath, dataname, frac, scale )

% impath = the root dir of the images

if ~exist( 'scale', 'var' ) || isempty( scale )
    scale = 1.0;
end

files = dir( [impath, '\*.bmp'] );
N = length( files );
img = imresize( imread( [impath, '\', files(1).name], 'bmp' ), scale );
[ n, m ] = size( img );
X = zeros( n, m, N );

for i = 1:N
    img = imresize( imread( [impath, '\', files(i).name], 'bmp' ), scale );
    X( :, :, i ) = double(img) / 255;
end

X = tensor( X );
[ b, Omega, linInd ] = sample_indices( X, frac );

data.X = X;
data.b = b;
data.Omega = Omega;
data.linInd = linInd;

save( ['..\..\data\',dataname], 'data' );
end