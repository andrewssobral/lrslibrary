function disp(a,name)
%DISP Command window display of a sptenmat.
%
%   DISP(T) displays the tensor without printing its name.
%
%   DISP(T,NAME) displays the tensor with the given name.
%
%   See also SPTENMAT/DISPLAY.
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


if ~exist('name','var')
    name = 'ans';
end

% Extract the number of nonzeros and number of dimensions
nz = size(a.vals,1);

% Print an intro sentence giving the name and the size
if (nz == 0)
    fprintf('%s is an all-zero sptenmat from an sptensor of size %s\n',...
        name, tt_size2str(a.tsize));
else
    fprintf('%s is a sptenmat from an sptensor of size %s with %d nonzeros\n',...
        name, tt_size2str(a.tsize), nz);
end

fprintf(1,'\t%s.rindices = %s (modes of tensor corresponding to rows)\n',...
        name,['[ ' num2str(a.rdims) ' ]'] );
fprintf(1,'\t%s.cindices = %s (modes of tensor corresponding to columns)\n',...
        name,['[ ' num2str(a.cdims) ' ]'] );


% Stop insane printouts
if (nz > 1000)
    r = input('Are you sure you want to print all nonzeros? (Y/N) ','s');
    if upper(r) ~= 'Y', return, end;
end

% Return now if there are no nonzeros
if (nz == 0)
    return;
end

% Pre-allocate memory for the output
output = cell(nz,1);
spc = floor(log10(max(a.subs)))+1;
fmt = ['\t(%' num2str(spc(1)) 'd,%' num2str(spc(2)) 'd)\t%g'];

for i = 1:nz
    output{i} = sprintf(fmt, a.subs(i,:), a.vals(i));
end
fprintf('%s\n',output{:});



