function z = isequal(x,y)
%ISEQUAL for sptensors.
%
%   ISEQUAL(A,B) compares the sparse tensors A and B for equality.
%
%   See also SPTENSOR.
%
%MATLAB Tensor Toolbox.
%Copyright 2012, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.
% http://www.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2012) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in the file LICENSE.txt


%% Observations for sparse matrix case.
% The result of isequal(a,full(a)) is true!

%%
if ~isequal(x.size,y.size) 
    z = false;
elseif isa(x,'sptensor') && isa(y,'sptensor') 
    z = (nnz(x-y) == 0);
elseif isa(y,'tensor')
    z = isequal(full(x),y);
else
    z = false;
end
