%load hall1-200; 
%[L,S,G,error,time]=GreBackground(XO,2,8,3,[144,176],1e-3,1);

%  load shoppingmall1-200; 
%  [L,S,G,error,time]=GreBackground(XO,3,6,3,[256,320],1e-3,1);
% 
load bootstrap1-200; 
[L,S,G,error,time]=GreBackground(XO,3,7,5,[120,160],1e-3,1);
% 
% load lobby1-200; 
% [L,S,G,error,time]=GreBackground(XO,2,6,3,[128,160],1e-3,1);