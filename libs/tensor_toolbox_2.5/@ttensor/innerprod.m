function res = innerprod(X,Y)
%INNERPROD Efficient inner product with a ttensor.
%
%   R = INNERPROD(X,Y) efficiently computes the inner product between
%   two tensors X and Y.  If Y is a tensor or sptensor, the inner
%   product is computed directly and the computational complexity
%   is...  If Y is a ktensor, the inner product method for that type
%   of tensor is called.
%
%   See also TTENSOR, KTENSOR/INNERPROD
%
%MATLAB Tensor Toolbox.
%Copyright 2010, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by Brett Bader and Tamara Kolda.
% http://csmr.ca.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2010) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in tensor_toolbox/LICENSE.txt
% $Id: innerprod.m,v 1.13 2010/03/19 23:46:31 tgkolda Exp $

% X is a ttensor
switch class(Y)
    
    case {'ttensor'}
        if ~isequal(size(X),size(Y))
            error('X and Y must be the same size.');
        end
        if prod(size(X.core)) > prod(size(Y.core))
            % Reverse argument and call this function again so that the
            % tensor with the smaller core is the first argument. 
            res = innerprod(Y,X);
            return
        end
        W = cell(ndims(X),1);
        for n = 1:ndims(X)
            W{n} = X.u{n}'*Y.u{n};
        end
        J = ttm(Y.core, W);
        res = innerprod(X.core,J);
        return

    case {'tensor','sptensor'}
        if ~isequal(size(X),size(Y))
            error('X and Y must be the same size.');
        end
        if (prod(size(X)) < prod(size(X.core)))
            Z = full(X);
            res = innerprod(Z,Y);
            return;
        end
        Z = ttm(Y,X.u,'t');
        res = innerprod(Z, X.core);
        return

    case {'ktensor'}
        % Reverse arguments to call ktensor implementation
        res = innerprod(Y,X);
        return
        
    otherwise
        disp(['Inner product not available for class ' class(Y)]);
        return
end


