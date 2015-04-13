function b = multitransp(a, dim)
%MULTITRANSP  Transposing arrays of matrices.
%    B = MULTITRANSP(A) is equivalent to B = MULTITRANSP(A, DIM), where
%    DIM = 1.
%
%    B = MULTITRANSP(A, DIM) is equivalent to
%    B = PERMUTE(A, [1:DIM-1, DIM+1, DIM, DIM+2:NDIMS(A)]), where A is an
%    array containing N P-by-Q matrices along its dimensions DIM and DIM+1,
%    and B is an array containing the Q-by-P transpose (.') of those N
%    matrices along the same dimensions. N = NUMEL(A) / (P*Q), i.e. N is
%    equal to the number of elements in A divided by the number of elements
%    in each matrix.
%
%    MULTITRANSP, PERMUTE and IPERMUTE are a generalization of TRANSPOSE
%    (.') for N-D arrays.
%
%    Example:
%       A 5-by-9-by-3-by-2 array may be considered to be a block array
%       containing ten 9-by-3 matrices along dimensions 2 and 3. In this
%       case, its size is so indicated:  5-by-(9-by-3)-by-2 or 5x(9x3)x2.
%       If A is ................ a 5x(9x3)x2 array of 9x3 matrices,
%       C = MULTITRANSP(A, 2) is a 5x(3x9)x2 array of 3x9 matrices.
%
%    See also PERMUTE, IPERMUTE, MULTIPROD, MULTITRACE, MULTISCALE.

% $ Version: 1.0 $
% CODE      by:                 Paolo de Leva (IUSM, Rome, IT) 2005 Sep 9
% COMMENTS  by:                 Code author                    2006 Nov 21
% OUTPUT    tested by:          Code author                    2005 Sep 13
% -------------------------------------------------------------------------

% Setting DIM if not supplied.
if nargin == 1, dim = 1; end

% Transposing
order = [1:dim-1, dim+1, dim, dim+2:ndims(a)];
b = permute(a, order);
