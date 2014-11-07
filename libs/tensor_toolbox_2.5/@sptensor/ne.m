function z = ne(x,y)
%NE Not equal (~=) for sptensors.
%
%   A ~= B compares the elements of A and B for equality. The arguments can
%   be a pair of sptensors, an sptensor and a tensor, or an sptensor and a
%   scalar.  Regardless, the result is always returned as a sparse tensor.
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
% The result of a ~= 5 is sparse.
% The result of a ~= 0 is sparse.
% The result of a ~= full(a) is sparse.

%% Case 1: One argument is a scalar
if isscalar(y)   
    if y == 0        
        z = sptensor(x.subs,true,size(x));
    else
        subs1 = x.subs(x.vals ~= y,:);
        subs2 = setdiff(allsubs(x),x.subs,'rows');
        z = sptensor([subs1;subs2],true,size(x));
    end   
    return;
end

% Call back with the arguments reversed.
if isscalar(x)
    z = ne(y,x);
    return;
end

%% Case 2: Both x and y are tensors or some sort
% Check that the sizes match
if ~isequal(x.size,y.size)
    error('Size mismatch');
end

% Case 2a: Two sparse tensors
if isa(x,'sptensor') && isa(y,'sptensor')

    % find entries where either x *or* y is nonzero, but not both
    subs1 = setxor(x.subs,y.subs,'rows'); 
    % find entries where both are nonzero, but inequal
    subs2 = intersect(x.subs,y.subs,'rows');
    subs2 = subs2(extract(x,subs2) ~= extract(y,subs2),:);
    % put it all together
    z = sptensor([subs1;subs2],true,size(x));
    return;

end

% Case 2b: y is a dense tensor
if isa(y,'tensor')
   
    % find entries where x is zero but y is nonzero
    subs1 = setdiff(allsubs(x),union(x.subs,find(y == 0),'rows'),'rows');
    
    % find entries where x is nonzero but not equal to y
    subs2 = x.subs(x.vals ~= y(x.subs,'extract'),:);

    % put it all together
    z = sptensor([subs1;subs2],true,size(x));
    return;
    
end


%% Otherwise
error('The arguments must be two sptensors or an sptensor and a scalar.');
