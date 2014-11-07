function data = tenInd2matInd_core( data, mode )

tcube = tenzeros( size(data.T) );
tcube( data.linInd ) = 1;
tmat = tenmat( tcube, mode );
data.matInd = find( double(tmat) );
data.rpca_mode = mode;

end