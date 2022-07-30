function S = mysign(A)
%MYSIGN True sign function with MYSIGN(0) = 1.

%   Called by various matrices in elmat/private.
%
%   Nicholas J. Higham
%   Copyright 1984-2005 The MathWorks, Inc.
%   $Revision: 1.4.4.1 $  $Date: 2005/11/18 14:15:15 $

S = sign(A);
S(find(S==0)) = 1;
