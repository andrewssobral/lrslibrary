function [v]=nsetdiff(A,B);
%NSETDIFF
%
%[v]=nsetdiff(A,B);
%Slow setdiff by CA, 1998

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

AS=sort(A);
len_AS=length(AS);
for i=1:len_AS-1,
    for j=i+1:len_AS,
        if AS(i)==AS(j),
            AS(i)=NaN;
        end;
    end
end;
I=find(isnan(AS));
if ~isempty(I)
    AS(I)=[];
end;


BS=sort(B);
len_BS=length(BS);
for i=1:len_BS-1,
    for j=i+1:len_BS,
        if BS(i)==BS(j),
            BS(i)=NaN;
        end;
    end
end;
I=find(isnan(BS));
if ~isempty(I)
    BS(I)=[];
end;


len_AS=length(AS);
len_BS=length(BS);
if len_AS >= len_BS
    for i=1:len_AS,
        for j=1:len_BS,
            if AS(i)==BS(j),
                AS(i)=NaN;
            end;
        end;
    end;
    I=find(isnan(AS));
    if ~isempty(I)
        AS(I)=[];
    end;
    v=AS;
else
    for i=1:len_BS,
        for j=1:len_AS,
            if BS(i)==AS(j),
                BS(i)=NaN;
            end;
        end;
    end;
    I=find(isnan(BS));
    if ~isempty(I)
        BS(I)=[];
    end;
    v=BS;
end;



