function export_data(A, fname)
%EXPORT_DATA Export tensor-related data to a file.
%   
%   EXPORT(A,FNAME) exports object A to the file named FNAME in plain ASCII
%   text. Export currently supports exporting the following data types:
%
%      - tensor
%      - matrix
%
%   In the case of a tensor, the first three lines give details about the
%   tensor. The format for a 4 x 3 x 2 tensor is as follows...
%
%      tensor
%      3
%      4 3 2
%      <value of A(1,1,1)>
%      <value of A(2,1,1)>
%      <value of A(3,1,1)>
%      <value of A(4,1,1)>
%      <value of A(1,2,1)>
%      <value of A(2,2,1)>
%      <value of A(3,2,1)>
%      <value of A(4,2,1)>
%      ...
%      <value of A(2,3,2)>
%      <value of A(3,3,2)>
%      <value of A(4,3,2)>
%
%   A matrix is formatted the same as a 2-way tensor except that the first
%   line says "matrix" rather than "tensor".
%
%   See also TENSOR, IMPORT_DATA
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

%% Open file
fid = fopen(fname,'w');
if (fid == -1)
    error('Cannot open file %s',fname);
end

%% Export the object

if isa(A,'tensor')
    
    fprintf(fid, 'tensor\n');
    export_size(fid, size(A));
    export_array(fid, A.data);   
 
elseif isnumeric(A) && ndims(A) == 2        

    fprintf(fid, 'matrix\n');
    export_size(fid, size(A));
    export_array(fid, A);    
    
else   
    
    error('Invalid data type for export');    
    
end


%% Close file
fclose(fid);

function export_size(fid, sz)
% Export the size of something to a file
fprintf(fid, '%d \n', length(sz)); % # of dimensions on one line
fprintf(fid, '%d ', sz); % # size of each dimensions on the next line
fprintf(fid, '\n');

function export_array(fid, data)
% Export dense data that supports numel and linear indexing
for i = 1:numel(data)
    fprintf(fid, '%.16e\n', data(i));
end
