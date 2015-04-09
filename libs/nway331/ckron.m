function C=ckron(A,B)
%CKRON
% C=ckron(A,B)
%
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

% Claus Andersson, Jan. 1996

% Should not be compiled to overwrite ckron.mex

[mA,nA] = size(A);
[mB,nB] = size(B);

C = zeros(mA*mB,nA*nB);
if mA*nA <= mB*nB
  for i = 1:mA
  iC = 1+(i-1)*mB:i*mB;
    for j = 1:nA
      jC = 1+(j-1)*nB:j*nB;
      C(iC,jC) = A(i,j)*B;
    end
  end
else
  for i = 1:mB
    iC = i:mB:(mA-1)*mB+i;
    for j = 1:nB
      jC = j:nB:(nA-1)*nB+j;
      C(iC,jC) = B(i,j)*A;
    end
  end
end
