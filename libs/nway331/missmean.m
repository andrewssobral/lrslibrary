function mm=missmean(X)

%MISSMEAN mean of a matrix X with NaN's
%
%[mm]=missmean(X)
%
%This function calculates the mean of a matrix X.
%X may hold missing elements denoted by NaN's which
%are ignored (weighted to zero).
%
%Check that for no column of X, all values are missing

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


%Insert zeros for missing, correct afterwards
missidx = isnan(X);
i = find(missidx);
X(i) = 0;

%Find the number of real(non-missing objects)
if min(size(X))==1,
   n_real=length(X)-sum(missidx);
else
   n_real=size(X,1)-sum(missidx);
end

i=find(n_real==0);
if isempty(i) %All values are real and can be corrected
   mm=sum(X)./n_real;
else %There are columns with all missing, insert missing
   n_real(i)=1;
   mm=sum(X)./n_real;
   mm(i)=i + NaN;
end
