function t = sptensor(varargin)
%SPTENSOR Create a sparse tensor.
%
%   X = SPTENSOR(SUBS, VALS, SZ, FUN) uses the rows of SUBS and VALS
%   to generate a sparse tensor X of size SZ = [m1 m2 ... mn]. SUBS is
%   an p x n array specifying the subscripts of the values to be
%   inserted into S. The k-th row of SUBS specifies the subscripts for
%   the k-th value in VALS. The values are accumulated at repeated
%   subscripts using the function FUN, which is specified by a
%   function handle.
%
%   There are several simplifications of this four argument call.
%
%   X = SPTENSOR(SUBS,VALS,SZ) uses FUN=@SUM.
%
%   X = SPTENSOR(SUBS,VALS) uses SM = max(SUBS,[],1).
%
%   X = SPTENSOR(SZ) abbreviates X = SPTENSOR([],[],SZ).
%
%   X = SPTENSOR(Y) copies/converts Y if it is an sptensor, an sptenmat, or
%   a dense tensor or MDA (the zeros are squeezed out), an sptensor3, or a
%   sparse matrix. Note that a row-vector, integer MDA is interpreted as a
%   size (see previous constructor).
%
%   S = SPTENSOR is the empty constructor.
%
%   The argument VALS may be scalar, which is expanded to be the
%   same length as SUBS, i.e., it is equivalent to VALS*(p,1).
%
%   Examples
%   subs = [1 1 1; 1 1 3; 2 2 2; 4 4 4; 1 1 1; 1 1 1]
%   vals = [0.5; 1.5; 2.5; 3.5; 4.5; 5.5]
%   siz = [4 4 4];
%   X = sptensor(subs,vals,siz) %<-- sparse 4x4x4, repeats summed
%   X = sptensor(subs,1,siz) %<-- scalar 2nd argument
%   X = sptensor(subs,vals,siz,@max) %<-- max for accumulation
%   myfun = @(x) sum(x) / 3;
%   X = sptensor(subs,vals,siz,myfun) %<-- custom accumulation
%
%   See also SPTENRAND, TENSOR, SPTENMAT, SPTENSOR3, ACCUMARRAY
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


% EMPTY Constructor
if (nargin == 0) || ((nargin == 1) && isempty(varargin{1}))
    t.subs = [];
    t.vals = [];
    t.size = [];
    t = class(t,'sptensor');
    return;
end

% SINGLE ARGUMENT
if (nargin == 1)

    source = varargin{1};

    switch(class(source))

        % COPY CONSTRUCTOR
        case 'sptensor',                
            t.subs = source.subs;
            t.vals = source.vals;
            t.size = source.size;
            t = class(t, 'sptensor');
            return;
            
        % CONVERT SPTENMAT
        case 'sptenmat',                

            % Extract the tensor size and order
            siz = source.tsize;

            % Convert the 2d-subscipts into nd-subscripts
            if ~isempty(source.rdims)
                subs(:,source.rdims) = ...
                    tt_ind2sub(siz(source.rdims),source.subs(:,1));
            end
            if ~isempty(source.cdims)
                subs(:,source.cdims) = ...
                    tt_ind2sub(siz(source.cdims),source.subs(:,2));
            end

            % Copy the values (which do not need to be modified)
            vals = source.vals;

            % Store everything
            t.subs = subs;
            t.vals = vals;
            t.size = siz;
            t = class(t, 'sptensor');
            return;

        % CONVERT TENSOR
        case 'tensor',                  
            [subs,vals] = find(source); 
            t.subs = subs;
            t.vals = vals;    
            t.size = size(source);
            t = class(t, 'sptensor');
            return;

        % CONVERT SPTENSOR3
        case 'sptensor3',
            K = size(source,3);
            [I,J] = size(source{1});
            nz = nnz(K);  
            bigsubs = [];
            bigvals = [];
            for k = 1:K
                [subs,vals] = find(source{k});
                if isempty(bigsubs)                 
                    bigsubs = [subs, k*ones(size(subs,1),1)];
                    bigvals = [vals];
                else
                    bigsubs = [bigsubs; subs, k*ones(size(subs,1),1)];
                    bigvals = [bigvals; vals];
                end
            end
            t.subs = bigsubs;
            t.vals = bigvals;
            t.size = [ I J K ];
            t = class(t,'sptensor');               
            return;
            
        % SPARSE MATRIX, SIZE, or MDA
        case {'numeric','logical','double'},

            % Case 1: SPARSE MATRIX
            if issparse(source)
                [i,j,s] = find(source);
                siz = size(source);
                t.subs = [i,j];
                t.vals = s;
                t.size = siz;
                t = class(t,'sptensor');
                return;
            end

            % Case 2: SPECIFYING THE SIZE
            if tt_sizecheck(source)
                t.subs = [];
                t.vals = [];
                t.size = source;
                t = class(t, 'sptensor');
                return;
            end

            % Case 3: An MDA
            t = sptensor(tensor(source));
            return;

    end % switch

end % nargin == 1

% SPECIAL CASE for INTERACTION WITH MEX FILES
if (nargin == 4) && (isnumeric(varargin{4})) && (varargin{4} == 0)

    % Store everything
    t.subs = varargin{1};
    t.vals = varargin{2};
    t.size = varargin{3};

    % Create the tensor
    t = class(t, 'sptensor');

    return;

end

% CONVERT A SET OF INPUTS
if (nargin == 2) || (nargin == 3) || (nargin == 4)

    % Extract the subscripts and values
    subs = varargin{1};   
    vals = varargin{2};
    
    tt_subscheck(subs);
    tt_valscheck(vals);
    if ~isempty(vals) && (numel(vals) ~= 1) && (size(vals,1) ~= size(subs,1))
        error('Number of subscripts and values must be equal');
    end

    % Extract the size
    if (nargin > 2)
        siz = varargin{3};
        tt_sizecheck(siz);
    else
        siz = max(subs,[],1);
    end

    % Check for wrong input
    if size(subs,2) > size(siz,2)
        error('More subscripts than specified by size')
    end

    % Check for subscripts out of range
    for j = 1:numel(siz)
        if ~isempty(subs) && max(subs(:,j)) > siz(j)
            error('Subscript exceeds sptensor size')
        end
    end

    % Extract the 'combiner' function handle
    if (nargin == 4)
        fun = varargin{4};
    else
        fun = @sum;
    end
    
    if isempty(subs)
        newsubs = [];
        newvals = [];
    else
        % Identify only the unique indices
        [newsubs,junk,loc] = unique(subs,'rows');

        % Sum the corresponding values
        newvals = accumarray(loc,vals,[size(newsubs,1) 1],fun);
    end

    % Find the nonzero indices of the new values
    nzidx = find(newvals);
    newsubs = newsubs(nzidx,:);
    newvals = newvals(nzidx);

    % Store everything
    t.subs = newsubs;
    t.vals = newvals;
    t.size = siz;

    % Create the tensor
    t = class(t, 'sptensor');

    return;
end

error('Unsupported use of function SPTENSOR.');

