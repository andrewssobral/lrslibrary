function [f,g] = tt_cp_fun(x,Z,Znormsqr)
%TT_CP_FUN Calculate function and gradient for CP fit function.
%
%  [F,G] = TT_CP_FUN(X,Z) where X is a vector containing the entries of the
%  components of the model and Z is the tensor to be fit.
%
%  See also TT_CP_VEC_TO_FAC, TT_FAC_TO_VEC, TT_CP_FG, CP_OPT
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


%% Convert x to a cell array of matrices
A = tt_cp_vec_to_fac(x,Z);

%% Call cp_fit and cp_gradient using cp_fg
[f,G] = tt_cp_fg(Z,A,Znormsqr);

%% Convert a cell array to a vector
g = tt_fac_to_vec(G);


