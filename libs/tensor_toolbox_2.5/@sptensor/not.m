function y = not(x)
%NOT Logical NOT (~) for sptensors.
%
%   ~X performs a logical not on the input tensor X. The result always
%   returned as a sparse tensor.
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
% The result of ~a is sparse.

%% Then compute those indicies that are not in x
subs = setdiff(allsubs(x),x.subs,'rows');

%% Assemble final result
y = sptensor(subs,true,x.size);
