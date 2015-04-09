function [DIA_Rel,DIA_Abs]=CoreDian(C)
%COREDIAN Calculates core diagonality
%
%[DIA_Rel,DIA_Abs]=CoreDian(C,Fac)
%
%DIA_Rel : Diagonality in percent
%DIA_Abs : Absolute sum of squares of
%          the diagonal elements 

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

Fac = size(C);
C = reshape(C,Fac(1),prod(Fac(2:end)));


c=size(Fac,2);
d=ones(1,c);
l=min(Fac);

C=C.^2;

DIA_Abs=0;
for j=1:l,
  [i,j]=getindxn(Fac,j*d);
  DIA_Abs = DIA_Abs + C(i,j);
end;

DIA_Rel=100*DIA_Abs/sum(sum(C));
