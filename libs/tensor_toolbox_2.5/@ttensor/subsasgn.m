function t = subsasgn(t,s,b)
%SUBSASGN Subscripted reference for a ttensor.
%
%   See also TTENSOR.
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
% $Id: subsasgn.m,v 1.9 2010/03/19 23:46:31 tgkolda Exp $

switch s(1).type
    case '.'
        switch s(1).subs
            case {'core','lambda'}
                if length(s) == 1
                    t = ttensor(b, t.u);
                else
                    tmpcore = subsasgn(t.core, s(2:end), b);
                    t = ttensor(tmpcore, t.u);
                end
            case {'u','U'}
                if length(s) == 1
                    t = ttensor(t.core, b);
                else
                    tmpu = subsasgn(t.u, s(2:end), b);
                    t = ttensor(t.core, tmpu);
                end
            otherwise
                error(['Cannot change field ', s.subs, ' directly.']);
        end
    case '()'
        error('Cannot change individual entries in ttensor.')
    case '{}'
        new_s(1).type = '.';
        new_s(1).subs = 'u';
        new_s(2:length(s)+1) = s;
        t = subsasgn(t, new_s, b);
    otherwise
        error('Invalid subsasgn.');
end


