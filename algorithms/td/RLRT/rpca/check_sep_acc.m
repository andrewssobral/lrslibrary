function [ acc, acc1 ] = check_sep_acc( E, truth )

img = double( E( :, :, 100 ) );
img( img > 0.1 ) = 1;
img( img <= 0.1 ) = 0;
acc1 = sum( img(:) == truth(:) ) / (size(img,1) * size(img,2));

% base acc (when img is completely black)
img = zeros( size(img) );
acc0 = sum( img(:) == truth(:) ) / (size(img,1) * size(img,2));

acc = (acc1 - acc0) / (1-acc0);

end