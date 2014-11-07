function X = tq_idwt( WX, wname )

global WT;
ndim = length(size(WX));

switch ndim
    case 2
        midr = size(WX,1)/2;      midc = size(WX,2)/2;
        cA = WX( 1:midr, 1:midc );
        cH = WX( 1:midr, midc+1:end );
        cV = WX( midr+1:end, 1:midc );
        cD = WX( midr+1:end, midc+1:end );
        X = idwt2( cA, cH, cV, cD, wname );
    
    case 3
        wt = WT;
        midr = size(WX,1)/2;    midc = size(WX,2)/2;    midk = size(WX,3)/2;
        for i = 1:2
            for j = 1:2
                for k = 1:2
                    wt.dec{i,j,k} = WX( (i-1)*midr+1 : i*midr, (j-1)*midc+1 : j*midc, (k-1)*midk+1 : k*midk );
                end
            end
        end
        X = idwt3( wt );
end
end