function [E]=GSM(V);
%GSM orthogonalization
%
% [E]=GSM(V);
% GS   Gram-Schmidt Method for orthogonalisation
%      An orthonormal basis spanning the columns of V is returned in E.
% 
%      This algorithm does not use pivoting or any other
%      stabilization scheme. For a completely safe orthogonalization
%      you should use 'ORTH()' though is may take triple the time.
%      'GSM()' is optimized for speed and requies only minimum storage
%      during iterations. No check of rank is performed on V!
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


[m n]=size(V);

%Allocate space for the basis
E=zeros(m,n);

%The first basis vector is taken directly from V
s=sqrt(sum(V(:,1).^2));
E(:,1)=V(:,1)/s;

%Find the other basis vectors as orthogonals to
%the already determined basis by projection
for k=2:n,
  f=V(:,k)-E(:,1:(k-1))*(E(:,1:(k-1))'*V(:,k));
  s=sqrt(sum(f.^2));
  if s<eps,
    E(:,k)=0*f;   %set to zeros
  else
    E(:,k)=f/s;   %normalize
  end;
end;
