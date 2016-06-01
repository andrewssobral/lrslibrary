clear;clc;

rng('default');
rng(100);

X = rand(300, 500);
omega = X<0.05;
%nnz(omega)

S = rand(500, 100);

XS1 = X* S;
Xom = X; Xom(~omega) = 0;
Xoc = X; Xoc( omega) = 0;
XS2 = Xom * S + Xoc * S;

if nnz(XS1 ~= XS2)
    disp('identical');
else
    disp('different');
end