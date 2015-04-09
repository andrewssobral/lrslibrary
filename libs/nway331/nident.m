function id = nident(J,order,mag);

% NIDENT make 'identity' multi-way array
%
% Make 'identity' array of order given by "order" and dimension J.
% id = nident(J,order);
% if extra input vector, mag, is given the j'th superdiagonal will be
% equal to mag(j)

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


% $ Version 2.00 $ May 2001 $ Changed to array notation $ RB $ Not compiled $
% $ Version 1.02 $ Date 28. July 1998 $ Not compiled $

if nargin<3
  mag=ones(J,1);
end

id=zeros(J^order,1);
for f=1:J
  idd=f;
  for i=2:order
    idd=idd+J^(i-1)*(f-1);
  end
  id(idd)=mag(f);
end
id  = reshape(id,ones(1,order)*J);
