function [Y,D]=normit(X)
%NORMIT normalize
%
%[Y,D]=normit(X);
%Normalizes X, such that X  = Y *D
%It also applies that    X' = D'*Y'
%Always give X column wise


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

[a b]=size(X);

Y=zeros(a,b);

SS=sqrt(sum(X.^2));
for i=1:b,
  Y(:,i)=X(:,i)./SS(i);
end;

if nargout==2,
  D=(Y'*Y)\(Y'*X);
end;

