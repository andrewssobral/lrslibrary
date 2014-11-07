function T = shrinkage_t( X, alpha )
% soft-thresholding on tensor X

Xv = double( X );
Tv = shrinkage_v( Xv, alpha );
T = tensor( Tv );
end