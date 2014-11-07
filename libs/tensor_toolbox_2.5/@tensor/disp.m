function disp(X,name)
%DISP Command window display of a tensor.
%
%   DISP(X) displays a tensor with no name.
%
%   DISP(X,NAME) displays a tensor with the given name.
%
%   See also TENSOR, TENSOR/DISPLAY.
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

fprintf(1,'%s is a tensor of size %s\n',name,tt_size2str(X.size));

if isempty(X.data)
    fprintf(1,'\t%s = []\n',name);
end

s = shiftdim(num2cell(X.data,1:2),2);

for i = 1:numel(s)
    fprintf('\t%s',name);
    if ndims(X) == 1
        fprintf('(:)');
    elseif ndims(X) == 2
        fprintf('(:,:)');
    elseif ndims(X) > 2
        fprintf('(:,:');
        fprintf(',%d',tt_ind2sub(X.size(3:end),i));
        fprintf(')');
    end
    fprintf(' = \n');
    output = tt_matrix2cellstr(s{i});
    fprintf('\t%s\n',output{:});
end

end
