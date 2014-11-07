function t = subsasgn(t,s,b)
%SUBSASGN Subscripted assignment for tenmat.  
%
%   Examples 
%   X = tenmat(rand(3,4,2),1); 
%   X(1:2,1:2) = ones(2,2); <-- Calls SUBSASGN 
%
%   See also TENMAT, TENMAT/SUBSREF.
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


switch s.type    
    case '()'
        [m n] = size(t.data);
	t.data(s.subs{:}) = b;
        if ~isequal([m n],size(t.data))
            error('Ambiguous change in size')
        end
    otherwise
        error('Invalid assignment for tenmat.')
end


