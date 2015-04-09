function [Y]=Refold3(X,R);
%REFOLD3
%
%[Y]=Refold3(X,R);
%
%If you have confused the way you unfolded
%your data, you can mend it by this function

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

Y=zeros(size(X));
Olis=zeros(1,R(3)*R(2));
lis=[0 kron([1:R(2)],R(3))];
lis=lis(1:R(2));
for k=1:R(3),
   Olis(1,(k-1)*R(2)+1:k*R(2))=k + lis;
end;
Y=X(1:R(1),Olis);