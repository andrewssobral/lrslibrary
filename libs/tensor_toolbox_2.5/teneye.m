function A = teneye(M,N)
%TENEYE Create identity tensor of specified size.
%
%   We say E is the "identity tensor" if TTSV(E,X,-1) = X for all X such
%   that NORM(X) = 1.
%
%   TENEYE(M,N) returns a (dense) identity tensor of order M and size N.
%   The identity tensor only exists when M is even and returns an error
%   otherwise. Due to the complexity of generating all possible indices,
%   this only works for relatively small M and N (i.e., N <= 7, M <= 6).
%
%   Examples
%   E = teneye(4,2); %<-- Create 2 x 2 x 2 x 2 identity tensor
%   x = rand(2,1); x = x/norm(x); %<-- Generate random x with norm 1
%   norm(ttsv(E,x,-1)-x) %<-- Check that ttsv(E,x,-1) = x
%
%   See also tensor, ttsv.
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

%% Check that it's an even tensor
if mod(M,2) ~= 0
    error('m must be even');
end

%% Generate all possible indices
idx = tt_combinator(N,M,'c','r');

%%
A = tenzeros(N*ones(1,M));
for i = 1:size(idx,1)
   p = perms(idx(i,:));
   for j = 1:M/2
       s(:,j) = (p(:,2*j-1) == p(:,2*j));
   end
   v = sum(sum(s,2)==M/2);
   A(p) = v / factorial(M);    
end

