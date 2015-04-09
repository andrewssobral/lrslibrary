function [Var_Rel,Var_Abs]=CoreVarn(C)
%COREVARN  Calculates the 'core-variance' 
%
%[Var_Rel,Var_Abs]=CoreVarn(C)
%Calculates the 'core-variance' using not means of
%variables(columns) but the avg of all elements in the matrix
%
%Var_Rel : Variance of squares in percent
%Var_Abs : Absolute variance of squares

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

C=C.^2;
p=prod(size(C));
mn=sum(C(:))/p;

Var_Abs=sum(sum( (C-mn).^2) );

Var_Rel=100*Var_Abs/(p*(p-1)*mn^2);


