function disp(t,name)
%DISP Command window display of a matricized tensor (tenmat).
%
%   DISP(T) displays a tensor as matrix with no name.
%
%   DISP(T,NAME) display a tensor as matrix with the given name.
%
%   See also TENMAT, TENMAT/DISPLAY.
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
    
fprintf('%s is a matrix corresponding to a tensor of size %s\n',...
        name,tt_size2str(t.tsize));
fprintf('\t%s.rindices = %s (modes of tensor corresponding to rows)\n',...
        name,['[ ' num2str(t.rindices) ' ]']);
fprintf('\t%s.cindices = %s (modes of tensor corresponding to columns)\n',...
        name,['[ ' num2str(t.cindices) ' ]']);

if isempty(t.data)
    fprintf('\t%s.data = []\n',name);
else
    fprintf('\t%s.data = \n',name);
    output = tt_matrix2cellstr(t.data);
    fprintf('\t\t%s\n',output{:});
end    


