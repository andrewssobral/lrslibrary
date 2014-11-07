function disp(t, name)
%DISP Command window display for a ktensor.
%
%   DISP(T) displays a Kruskal tensor with no name.
%
%   DISP(T,NAME) display a Kruskal tensor with the given name.
%
%   See also DISP, KTENSOR/DISPLAY, KTENSOR
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

fprintf('%s is a ktensor of size %s\n', name, tt_size2str(size(t)));
fprintf('\t%s.lambda = %s\n',name, ['[ ' num2str(t.lambda') ' ]'] );

if (ndims(t) > 0)
    for j = 1 : ndims(t)
        fprintf('\t%s.U{%d} = \n', name, j);
        output = tt_matrix2cellstr(t.u{j});
        fprintf('\t\t%s\n',output{:});
    end
end
