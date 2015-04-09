function [X]=missmult(A,B)

%MISSMULT product of two matrices containing NaNs
%
%[X]=missmult(A,B)
%This function determines the product of two matrices containing NaNs
%by finding X according to
%     X = A*B
%If there are columns in A or B that are pur missing values,
%then there will be entries in X that are missing too.
%
%The result is standardized, that is, corrected for the lower
%number of contributing terms.
%
%Missing elements should be denoted by 'NaN's


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

%INBOUNDS
%REALONLY

[ia ja]=size(A);
[ib jb]=size(B);
X=zeros(ia,jb);

one_arry=ones(ia,1);
for j=1:jb,
   p=one_arry*B(:,j)';
   tmpMat=A.*p;
   X(:,j)=misssum(tmpMat')';
end;
