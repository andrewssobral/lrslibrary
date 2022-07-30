function WX = tq_dwt( X, wname )
% The dimensions of X have to be even!!

ndim = length(size(X));

switch ndim
    case 2
        [ cA, cH, cV, cD ] = dwt2( X, wname );
        WX = [ cA, cH; cV, cD ];

    case 3
        WX = zeros( size(X) );
        midk = size(X,3) / 2;
        wt = dwt3( X, wname );
        WX(:,:,1:midk) = [ wt.dec{1,1,1}, wt.dec{1,2,1}; wt.dec{2,1,1}, wt.dec{2,2,1} ];
        WX(:,:,midk+1:end) = [ wt.dec{1,1,2}, wt.dec{1,2,2}; wt.dec{2,1,2}, wt.dec{2,2,2} ];
end

if sum( size(X) ~= size(WX) )
    disp( 'Dimensions of X have to be even!' );
end

end