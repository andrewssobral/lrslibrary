function data = tenInd2matInd( dataname, mode )
% convert the indices in Omega w.r.t. the data tensor to indices w.r.t. the
% mode-i unfolding of the tensor

datapath = ['..\..\data\',dataname];
load( datapath );
data = tenInd2matInd_core( data, mode );

save( datapath, 'data' );
end