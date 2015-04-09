function st=stdnan(X);

%STDNAN estimate std with NaN's
%
% Estimates the standard deviation of each column of X
% when there are NaN's in X.
%
% Columns with only NaN's get a standard deviation of zero


% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $

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


[I,J]=size(X);

st=[];
for j=1:J
  id=find(~isnan(X(:,j)));
  if length(id)
    st=[st std(X(id,j))];
  else
    st=[st 0];
  end
end