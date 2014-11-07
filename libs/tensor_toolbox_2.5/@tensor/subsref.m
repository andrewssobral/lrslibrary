function a = subsref(t,s)
%SUBSREF Subscripted reference for tensors.
%
%   We can extract elements or subtensors from a tensor in the
%   following ways.
%
%   Case 1a: y = X(i1,i2,...,iN), where each in is an index, returns a
%   scalar.
%
%   Case 1b: Y = X(R1,R2,...,RN), where one or more Rn is a range and
%   the rest are indices, returns a sparse tensor.
%
%   Case 2a: V = X(S) or V = X(S,'extract'), where S is a p x n array
%   of subscripts, returns a vector of p values.
%
%   Case 2b: V = X(I) or V = X(I,'extract'), where I is a set of p
%   linear indices, returns a vector of p values.
%
%   Any ambiguity results in executing the first valid case. This
%   is particularly an issue if ndims(X)==1.
%
%   Examples
%   X = tensor(rand(3,4,2,1),[3 4 2 1]);
%   X(1,1,1,1) %<-- produces a scalar
%   X(1,1,1,:) %<-- produces a tensor of order 1 and size 1
%   X(:,1,1,:) %<-- produces a tensor of size 3 x 1
%   X(1:2,[2 4],1,:) %<-- produces a tensor of size 2 x 2 x 1
%   X(1:2,[2 4],1,1) %<-- produces a tensor of size 2 x 2
%   X([1,1,1,1;3,4,2,1]) %<-- returns a vector of length 2
%   X = tensor(rand(10,1),10);
%   X([1:6]') %<-- extracts a subtensor
%   X([1:6]','extract') %<-- extracts a vector of 6 elements
%
%   See also TENSOR, TENSOR/FIND.
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


switch s(1).type
    case '{}'
        error('Cell contents reference from a non-cell array object.')
    case '.'
        fieldname = s(1).subs;
        switch fieldname
            case 'data'
                a = tt_subsubsref(t.data,s);
            case 'size'
                a = tt_subsubsref(t.size,s);
            otherwise
                error(['No such field: ', fieldname]);
        end
    case '()'

        % *** CASE 1: Rectangular Subtensor ***
        if (numel(s(1).subs) == ndims(t)) && ~isequal(s(1).subs{end},'extract')

            % Copy the subscripts
            region = s(1).subs;

            if numel(region) ~= ndims(t)
                error('Invalid number of subscripts');
            end

            % Extract the data
            newdata = t.data(region{:});

            % Determine the subscripts
            newsiz = [];                % (future) new size
            kpdims = [];                % dimensions to keep
            rmdims = [];                % dimensions to remove

            % Determine the new size and what dimensions to keep
            for i = 1:length(region)
                if ischar(region{i}) && (region{i} == ':')
                    newsiz = [newsiz size(t,i)];
                    kpdims = [kpdims i];
                elseif numel(region{i}) > 1
                    newsiz = [newsiz numel(region{i})];
                    kpdims = [kpdims i];
                else
                    rmdims = [rmdims i];
                end
            end

            % If the size is zero, then the result is returned as a scalar;
            % otherwise, we convert the result to a tensor.
            if isempty(newsiz)
                a = newdata;
            else
                if isempty(rmdims)
                    a = tensor(newdata,newsiz);
                else
                    a = tensor(permute(newdata,[kpdims rmdims]),newsiz);
                end
            end
            a = tt_subsubsref(a,s);
            return;
        end

        % *** CASE 2a: Subscript indexing
        if size(s(1).subs{1},2) == ndims(t)
            % extract array of subscripts
            subs = s(1).subs{1};
            a = t.data(tt_sub2ind(t.size,subs));
            a = tt_subsubsref(a,s);
            return;
        end

        % *** CASE 2b: Linear indexing ***
        if numel(s(1).subs) ~= 1
            error('Invalid indexing');
        end

        idx = s(1).subs{1};
        if ndims(idx) ~=2 || size(idx,2) ~= 1
            error('Expecting a column index');
        end

        a = t.data(idx);
        a = tt_subsubsref(a,s);
        return;
end
