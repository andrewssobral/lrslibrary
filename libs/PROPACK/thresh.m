function A = thresh(B,tol)
% THRESH converts all entries less than TOL equal to zero
%   A = THRESH(B) uses the default tolerance: 100*eps
%       (where eps is matlab's machine precision)
%   A = THRESH(B,TOL) uses the tolerance TOL
%
%  written by Stephen Becker, 11/17/07 srbecker@caltech.edu
%   See also eps

if nargin < 2  tol = 100*eps; end
if isreal(B) 
    A = B;
    A( find( abs(B) < tol ) ) = 0;
else
    A = thresh(real(B),tol) + i*thresh(imag(B),tol);
end