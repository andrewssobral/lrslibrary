function G=neye(Fac);
% NEYE  Produces a super-diagonal array
%
%function G=neye(Fac);
%
% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 1.00 $ Date 5. Aug. 1998 $ Not compiled $
%
% This algorithm requires access to:
% 'getindxn'
%
% See also:
% 'parafac' 'maxvar3' 'maxdia3'
%
% ---------------------------------------------------------
%             Produces a super-diagonal array
% ---------------------------------------------------------
%	
% G=neye(Fac);
%
% Fac      : A row-vector describing the number of factors
%            in each of the N modes. Fac must be a 1-by-N vector. 
%            Ex. [3 3 3] or [2 2 2 2]

% Copyright (C) 1995-2006  Rasmus Bro & Claus Andersson
% Copenhagen University, DK-1958 Frederiksberg, Denmark, rb@life.ku.dk
%
% This program is free software; you can redistribute it and/or modify it under 
% the terms of the GNU General Public License as published by the Free Software 
% Foundation; either version 2 of the License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT 
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS 
% FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
% You should have received a copy of the GNU General Public License along with 
% this program; if not, write to the Free Software Foundation, Inc., 51 Franklin 
% Street, Fifth Floor, Boston, MA  02110-1301, USA.

N=size(Fac,2);
if N==1,
   fprintf('Specify ''Fac'' as e vector to define the order of the core, e.g.,.\n')
   fprintf('G=eyecore([2 2 2 2])\n')
end;

G=zeros(Fac(1),prod(Fac(2:N)));

for i=1:Fac(1),
   [gi,gj]=getindxn(Fac,ones(1,N)*i);
   G(gi,gj)=1;
end;

G = reshape(G,Fac);