function [tf, tf_lambda, tf_U] = isequal(A,B)
%ISEQUAL True if each component of two ktensor's is numerically equal.
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

tf = false;
tf_lambda = false;
tf_U = false;

if ~isa(B,'ktensor')
    return;
end    

tf_lambda = isequal(A.lambda, B.lambda);
if ncomponents(A) == ncomponents(B)
    tf_U = cellfun(@isequal, A.u, B.u);
end
tf = tf_lambda & all(tf_U);

