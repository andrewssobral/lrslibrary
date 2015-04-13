function y = XonOmega(U,V,varargin)
% y = XonOmega(U,V,omega)
%   This implicitly forms the matrix A = U*V'
%   and then returns A(omega).
%
% y = XonOmega(U,V,I,J)
%   does the same thing, but uses subscript indices I and J
%   ( where [I,J] = ind2sub( size(A), omega )
%
% y = XonOmega(U,V,OMEGA)
%   does the same thing, but now omega is specified by the nonzero
%   entries of the matrix OMEGA.   This will give same results
%   as previous variations IF the vector omega is sorted.
%
%   Stephen Becker, 11/10/08


% This file only runs if there is no compiled mex file in the
% path.  You do NOT want to run this file as-is, because it's
% very wasteful, except when length(omega) is comparable to nnz(U*V').

str = ['Using slow matlab code for XonOmega because mex file not compiled\n',...
    'To disable this warning in the future, run the command:\n',...
    '   warning(''off'',''SVT:NotUsingMex'')   '];

warning('SVT:NotUsingMex',str);

M = size(U,1);
N = size(V,1);
FULLMATRIX =  ( M*N < 50*50  );
if nargin > 3 
    I = varargin{1};
    J = varargin{2};
    y = zeros(size(I));
    if length(I) == M*N, FULLMATRIX = true; end
    if FULLMATRIX
        omega = sub2ind( [M,N], I, J );
    end
else
    omega = varargin{1};
    if length(omega) == M*N, FULLMATRIX = true; end
    if isvector(omega)
        if ~FULLMATRIX
            [I,J] = ind2sub(  [M,N], omega );
            y = zeros(size(omega));
        end
    else
        if FULLMATRIX
            omega = find( omega );
            y = zeros(size(omega));
        else
            [I,J] = find( omega );
            y = zeros(size(I));
        end
    end
end


if FULLMATRIX
    % if Omega is the entire index set, or
    % if U and V are small, this is simplest and fastest:

    A = U*V';
    y = A(omega);
    disp('FULL MATRIX');

else
    % this will be slow, but for large matrices it's the only
    % practical option
    
    for k = 1:length(I)
        y(k) = U( I(k), :) * ( V( J(k), : )' );
    end
end
