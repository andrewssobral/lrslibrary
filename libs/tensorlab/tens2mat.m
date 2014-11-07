function M = tens2mat(T,mode_row,mode_col)
%TENS2MAT Matricize a tensor.
%   M = tens2mat(T,mode_row,mode_col) matricizes a tensor T into a matrix M
%   of dimensions prod(size_tens(mode_row))-by-prod(size_tens(mode_col)),
%   where size_tens is equal to size(T). The columns (rows) of M are
%   obtained by fixing the indices of T corresponding to mode_col
%   (mode_row) and looping over the remaining indices in the order mode_row
%   (mode_col). E.g., if A and B are two matrices and T = cat(3,A,B), then
%   tens2mat(T,1:2,3) is the matrix [A(:) B(:)].
%
%   M = tens2mat(T,mode_row) matricizes a tensor T, where mode_col is
%   chosen as the sequence [1:ndims(T)]\mode_row.
%
%   M = tens2mat(T,[],mode_col) matricizes a tensor T, where mode_row is
%   chosen as the sequence [1:ndims(T)]\mode_col.
%
%   See also mat2tens.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check arguments.
size_tens = size(T);
N = length(size_tens);
size_tens = [size_tens 1];
if nargin <= 2, mode_col = []; end
if isempty(mode_row) && isempty(mode_col)
    error('tens2mat:InvalidModes', ...
          'Either mode_row or mode_col must be non-empty.');
end
mode_row = mode_row(mode_row <= N); % >N is treated as singleton dimension.
mode_col = mode_col(mode_col <= N);
if isempty(mode_col), mode_col = complement(mode_row,N); end
if isempty(mode_row), mode_row = complement(mode_col,N); end
if isempty(mode_col), mode_col = N+1; end
if isempty(mode_row), mode_row = N+1; end

% Matricize the tensor.
if any(mode_row(:).' ~= 1:length(mode_row)) || ...
   any(mode_col(:).' ~= length(mode_row)+(1:length(mode_col)))
    T = permute(T,[mode_row(:).' mode_col(:).']);
end
M = reshape(T,prod(size_tens(mode_row)),[]);

end

function mode_col = complement(mode_row,N)
    bits = ones(1,N);
    bits(mode_row) = 0;
    mode_col = 1:N;
    mode_col = mode_col(logical(bits(mode_col)));
end
