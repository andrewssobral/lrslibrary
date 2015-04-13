function Y = multiskew(X)
% Returns the skew-symmetric parts of the matrices in the 3D matrix X.
%
% function Y = multiskew(X)
%
% Y is a 3D matrix the same size as X. Each slice Y(:, :, i) is the
% skew-symmetric part of the slice X(:, :, i).
%
% See also: multiprod multitransp multiscale multisym

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Jan. 31, 2013.
% Contributors: 
% Change log: 

    Y = .5*(X - multitransp(X));
    
end
