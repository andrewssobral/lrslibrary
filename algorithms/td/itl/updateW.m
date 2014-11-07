%update w(t) to w(t+1)
%para:
%j: index of output y
%old_w: W matrix
%X: the input vector X(t) 
%lambda: forgetting factor between [0,1]
%----------------------------------------------------------------
% Copyright: 2005,
%         Spiros Papadimitriou, Jimeng Sun, Christos Faloutsos.
%         All rights reserved.
% Please address questions or comments to: jimeng@cs.cmu.edu
%----------------------------------------------------------------
function [w,d,x]=updateW(old_x, old_w,old_d, lambda)
y=old_w'*old_x;
d=lambda*old_d+y^2;
e=old_x-old_w*y;
w=old_w+e*y/d;
x=old_x-w*y;
w=w/norm(w);