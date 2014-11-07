function A = import_data(fname)
%IMPORT_DATA Import tensor-related data to a file.
%
%   A = IMPORT_DATA(FNAME) imports an object A from the file named FNAME.
%   The supported data types and formatting of the file are explained in
%   EXPORT_DATA. 
%
%   See also TENSOR, EXPORT_DATA
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
fid = fopen(fname,'r');
if (fid == -1)
    error('Cannot open file %s',fname);
end

%% Get the type of object
line = fgets(fid);
type = sscanf(line, '%s');

%% Import the object

if strcmpi(type,'tensor')
    
    sz = import_size(fid);
    data = import_array(fid, prod(sz));   
    A = tensor(data, sz);
 
elseif strcmpi(type,'matrix')        

    sz = import_size(fid);
    data = import_array(fid, prod(sz));   
    A = reshape(data, sz);
    
else   
    
    error('Invalid data type for export');    
    
end


%% Close file
fclose(fid);

function sz = import_size(fid)
% Import the size of something from a file
line = fgets(fid);
n = sscanf(line, '%d');
line = fgets(fid);
sz = sscanf(line, '%d');
sz = sz';
if (size(sz,2) ~= n)
    error('Imported dimensions are not of expected size');
end

function data = import_array(fid, n)
% Export dense data that supports numel and linear indexing
data = fscanf(fid, '%e', n);
