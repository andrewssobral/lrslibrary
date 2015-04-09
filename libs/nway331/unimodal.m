function B=unimodal(X,Y,Bold)

%UNIMODAL unimodal regression
%
% Solves the problem min|Y-XB'| subject to the columns of 
% B are unimodal and nonnegative. The algorithm is iterative
% If an estimate of B (Bold) is given only one iteration is given, hence
% the solution is only improving not least squares
% If Bold is not given the least squares solution is estimated
%
% I/O B=unimodal(X,Y,Bold)
%
% Reference
% Bro and Sidiropoulos, "Journal of Chemometrics", 1998, 12, 223-247. 



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


if nargin==3
   B=Bold;
   F=size(B,2);
   for f=1:F
     y=Y-X(:,[1:f-1 f+1:F])*B(:,[1:f-1 f+1:F])';
     beta=pinv(X(:,f))*y;
     B(:,f)=ulsr(beta',1);
   end
else
   F=size(X,2);
   maxit=100;
   B=randn(size(Y,2),F);
   Bold=2*B;
   it=0;
   while norm(Bold-B)/norm(B)>1e-5&it<maxit
     Bold=B;
     it=it+1;
     for f=1:F
       y=Y-X(:,[1:f-1 f+1:F])*B(:,[1:f-1 f+1:F])';
       beta=pinv(X(:,f))*y;
       B(:,f)=ulsr(beta',1);
     end
   end
   if it==maxit
     disp([' UNIMODAL did not converge in ',num2str(maxit),' iterations']);
   end
end