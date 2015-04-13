function Y = updateSparse_slow(Y,b,indx,i,j)
% This is a replacement for the mex-file that updates the values of a sparse matrix Y.
%   The mex-file should be called as:
%
%   updateSparse(Y,b)
%
%   which will implicitly do the following:  Y(omega) = b
%   where "omega" is the set of nonzero indices of Y (in linear
%   ordering, i.e. column-major ordering).
%
%   If "omega" is not sorted, then you must do the following:
%
%       [temp,indx] = sort(omega);  % we don't care about "temp"
%       updateSparse(Y,b,indx);
%
%   which will ensure that everything is in the proper order.
%
%   If you see a warning about missing mex files, it means the mex
%   files have not been pre-compiled for your system.
%   In the [SVT root]/private directory, run the "install_mex.m" script.
%
%   If that fails, then this program will execute instead.

% This file and mex file by Stephen Becker, srbecker@caltech.edu 11.12.08
% Modified 4.20.11 by Farshad Harirchi and Stephen Becker

str = ['Using slow matlab code for updateSparse because mex file not compiled\n',...
    'To disable this warning in the future, run the command:\n',...
    '   warning(''off'',''SVT:NotUsingMex'')\n ',...
    'or install the mex file by running:\n',...
    '   mex updateSparse.c'];

warning('SVT:NotUsingMex',str);

if nargin<5,
    [i,j,s] = find(Y);
end

if nargin > 2 && ~isempty(indx)
    % resort
    b = b(indx);
end

% [i,j,s] = find(Y);
[m,n] = size(Y);
Y = sparse(i,j,b,m,n);
