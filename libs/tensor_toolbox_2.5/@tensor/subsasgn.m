function x = subsasgn(x,s,b)
%SUBSASGN Subscripted assignment for a tensor.
%
%   We can assign elements to a tensor in three ways.
%
%   Case 1: X(R1,R2,...,RN) = Y, in which case we replace the
%   rectangular subtensor (or single element) specified by the ranges
%   R1,...,RN with Y. The right-hand-side can be a scalar, a tensor, or an
%   MDA. 
%
%   Case 2a: X(S) = V, where S is a p x n array of subscripts and V is
%   a scalar or a vector containing p values.
%
%   Case 2b: X(I) = V, where I is a set of p linear indices and V is a
%   scalar or a vector containing p values. Resize is not allowed in this
%   case. 
%
%   Examples
%   X = tensor(rand(3,4,2))
%   X(1:2,1:2,1) = ones(2,2) %<-- replaces subtensor
%   X([1 1 1;1 1 2]) = [5;7] %<-- replaces two elements
%   X([1;13]) = [5;7] %<-- does the same thing
%   X(1,1,2:3) = 1 %<-- grows tensor
%   X(1,1,4) = 1 %<- grows the size of the tensor
%
%   See also TENSOR, TENSOR/SUBSREF.
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


switch s.type
    case '.'
        error(['Cannot change field ', s.subs, ' directly.']);

    case '{}'
        error('Subscript cell reference not supported for tensor.');

    case '()'

        % We don't allow sub-subscriptingfor tensors.
        if numel(s) ~= 1
            error('Invalid subscripting');
        end

        % Figure out if we are doing a subtensor, a list of subscripts, or
        % a list of linear indices...
        type = 'error';
        if ndims(x) <= 1
           if (numel(s.subs) > 1) || ...
                   ((ndims(s.subs{1}) == 2) && ...
                    (size(s.subs{1},1) == 1  || size(s.subs{1},2) == 1))
               type = 'subtensor';
           elseif ndims(s.subs{1}) == 2
               type = 'subscripts';
           end
        else
            if numel(s.subs) >= ndims(x)
                type = 'subtensor';
            elseif ndims(s.subs{1}) == 2
                if size(s.subs{1},2) >= ndims(x)
                    type = 'subscripts';
                elseif size(s.subs{1},2) == 1
                    type = 'linear indicies';
                end
            end
        end


        % *** CASE 1: Rectangular Subtensor ***
        if isequal(type,'subtensor')
            if isa(b,'tensor')
                x.data(s.subs{:},1) = b.data;
            else
                x.data(s.subs{:},1) = b;
            end
            % Check if the size has grown!
            % Can't vectorize this due to possible trailing 1's
            for i = 1:numel(x.size)
                x.size(i) = max(x.size(i),size(x.data,i));
            end
            % Check if order has grown
            for i = numel(x.size)+1:numel(s.subs)
                x.size(i) = size(x.data,i);
            end
            return;
        end

        % *** CASE 2a: Subscript indexing ***
        if isequal(type,'subscripts');

            % extract array of subscripts
            subs = s.subs{1};

            % will the size change? if so, we first need to resize x
            n = ndims(x);
            bsiz = max(subs,[],1);
            newsiz = [max([x.size;bsiz(1:n)]) bsiz(n+1:end)];
            if ~isequal(newsiz,x.size)
                % We need to enlarge x.data. A trick is to assign its last
                % "new" element to zero. This resizes the array correctly.
                if numel(newsiz) == 1
                    str = sprintf('x.data(%d)=0;',newsiz);
                else
                    str = [sprintf('x.data(') ...
                        sprintf('%d,',newsiz(1:end-1)) ...
                        sprintf('%d)=0;', newsiz(end)) ];
                end
                eval(str);
                x.size = newsiz;
            end

            % finally, we can copy in the new data
            x.data(tt_sub2ind(newsiz,subs)) = b;
            return;
        end

        % *** CASE 2b: Linear indexing ***
        if isequal(type,'linear indicies');
            idx = s.subs{1};
            x.data (idx) = b;
            return
        end
        
        error('Invalid use of tensor/subsasgn');
end





