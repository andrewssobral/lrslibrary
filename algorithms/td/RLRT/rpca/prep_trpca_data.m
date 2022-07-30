function T = prep_trpca_data( impath, dataname, scale )

% impath = the root dir of the images

if ~exist( 'scale', 'var' ) || isempty( scale )
    scale = 1.0;
end

files = dir( [impath, '\*.bmp'] );
N = length( files );
% N = N - 1;
img = imresize( imread( [impath, '\', files(1).name], 'bmp' ), scale );
img = rgb2gray( img );
sz = size( img );   lsz = length(sz);
sz(lsz+1) = N;
X = zeros( sz );

for i = 1:N
    img = imresize( imread( [impath, '\', files(i).name], 'bmp' ), scale );
    img = rgb2gray( img );
    if lsz == 2
        X( :, :, i ) = double(img) / 255;
    else
        X( :, :, :, i ) = double(img) / 255;
    end
end

% img = imresize( imread( [impath, '\', files(N+1).name], 'bmp' ), scale );
% img = rgb2gray( img );
% data.truth = double(img) / 255;

T = tensor( X );

data.X = T;

save( ['..\..\data\',dataname], 'data' );
end