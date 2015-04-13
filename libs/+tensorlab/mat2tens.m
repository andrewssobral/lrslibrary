function T = mat2tens(M,size_tens,mode_row,mode_col)
%MAT2TENS Tensorize a matrix.
%   T = mat2tens(M,size_tens,mode_row,mode_col) tensorizes a matrix M into
%   a tensor T of dimensions size_tens, given its matricization defined by
%   mode_row and mode_col. The columns (rows) of M should correspond to
%   fixing the indices of T corresponding to mode_col (mode_row) and
%   looping over the remaining indices in the order mode_row (mode_col).
%   E.g., if A and B are two matrices and M = [A B], then
%   mat2tens(M,[size(A) 2],1,2:3) is the tensor T = cat(3,A,B).
%
%   T = mat2tens(M,size_tens,mode_row) tensorizes a matrix M, where
%   mode_col is chosen as the sequence [1:length(size_tens)]\mode_row.
%
%   T = mat2tens(M,size_tens,[],mode_col) tensorizes a matrix M, where
%   mode_row is chosen as the sequence [1:length(size_tens)]\mode_col.
%
%   See also tens2mat.

%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%            Marc Van Barel (Marc.VanBarel@cs.kuleuven.be)
%            Lieven De Lathauwer (Lieven.DeLathauwer@kuleuven-kulak.be)

% Check arguments.
N = length(size_tens);
size_tens = [size_tens 1];
if nargin <= 3, mode_col = []; end
if isempty(mode_row) && isempty(mode_col)
    error('mat2tens:InvalidModes', ...
          'Either mode_row or mode_col must be non-empty.');
end
mode_row = mode_row(mode_row <= N); % >N is treated as singleton dimension.
mode_col = mode_col(mode_col <= N);
if isempty(mode_col), mode_col = complement(mode_row,N); end
if isempty(mode_row), mode_row = complement(mode_col,N); end
if isempty(mode_col), mode_col = N+1; end
if isempty(mode_row), mode_row = N+1; end
if size(M,1) ~= prod(size_tens(mode_row)) || ...
   size(M,2) ~= prod(size_tens(mode_col))
    error('mat2tens:InvalidDimensions', ...
          'Invalid matrix dimensions for the requested tensorization.');
end

% Tensorize the matrix.
T = reshape(M,size_tens([mode_row(:).' mode_col(:).']));
if any(mode_row(:).' ~= 1:length(mode_row)) || ...
   any(mode_col(:).' ~= length(mode_row)+(1:length(mode_col)))
    iperm([mode_row(:).' mode_col(:).']) = ...
        1:length(mode_row)+length(mode_col);
    T = permute(T,iperm);
end

end

function mode_col = complement(mode_row,N)
    bits = ones(1,N);
    bits(mode_row) = 0;
    mode_col = 1:N;
    mode_col = mode_col(logical(bits(mode_col)));
end
