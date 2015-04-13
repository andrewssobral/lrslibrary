function tr = multitrace(A)
% Computes the traces of the 2D slices in a 3D matrix.
% 
% function tr = multitrace(A)
%
% For a 3-dimensional matrix A of size n-by-n-by-N, returns a column vector
% tr of length N such that tr(k) = trace(A(:, :, k));
%
% See also: multiprod multitransp multiscale

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Dec. 30, 2012.
% Contributors: 
% Change log: 

    
    assert(ndims(A) <= 3, ...
           ['multitrace is only well defined for matrix arrays of 3 ' ...
            'or less dimensions.']);

	tr = diagsum(A, 1, 2);

end
