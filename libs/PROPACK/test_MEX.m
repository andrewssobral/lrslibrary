%{
  Tests mex files

installing on windows, add this to compiler options:
rem SRB changing this: use /MT library instead of /MD for a static build
set COMPFLAGS=/c /Zp8 /GR /W3 /EHs /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0 /DMATLAB_MEX_FILE /nologo /MT
rem set COMPFLAGS=/c /Zp8 /GR /W3 /EHs /D_CRT_SECURE_NO_DEPRECATE /D_SCL_SECURE_NO_DEPRECATE /D_SECURE_SCL=0 /DMATLAB_MEX_FILE /nologo /MD

%}
disp('============== Testing mex files ===================');

disp('-- Testing XonOmega --');
M = 100; N = 60; R = 17;
U=randn(M,R); V = randn(N,R);
omega = randperm(M*N); omega = sort( omega(1:round(end/2) ) )';
A = U*V'; y0 = A(omega);
[I,J] = ind2sub(size(A), omega );
Y = sparse( I, J, ones(size(omega)) );
y1 = XonOmega(U,V,omega);   e1 = norm(y1-y0);
y2 = XonOmega(U,V,Y);       e2 = norm(y2-y0);
y3 = XonOmega(U,V,I,J);     e3 = norm(y3-y0);
fprintf('\tMethod 1 error: %e\n\tMethod 2 error: %e\n\tMethod 3 error: %e\n',e1,e2,e3);
disp('-- Testing XonOmegaTranspose --');
y1 = XonOmegaTranspose(U.',V.',omega);   e1 = norm(y1-y0);
y2 = XonOmegaTranspose(U.',V.',Y);       e2 = norm(y2-y0);
y3 = XonOmegaTranspose(U.',V.',I,J);     e3 = norm(y3-y0);
fprintf('\tMethod 1 error: %e\n\tMethod 2 error: %e\n\tMethod 3 error: %e\n',e1,e2,e3);


disp('-- Testing XonOmega with complex numbers--');
U=randn(M,R) + 1i*randn(M,R); V = randn(N,R) + 1i*randn(N,R);
A = U*V'; y0 = A(omega);
[I,J] = ind2sub(size(A), omega );
Y = sparse( I, J, ones(size(omega)) );
y1 = XonOmega(U,V,omega);   e1 = norm(y1-y0);
y2 = XonOmega(U,V,Y);       e2 = norm(y2-y0);
y3 = XonOmega(U,V,I,J);     e3 = norm(y3-y0);
fprintf('\tMethod 1 error: %e\n\tMethod 2 error: %e\n\tMethod 3 error: %e\n',e1,e2,e3);

disp('-- Testing XonOmegaTranspose with complex numbers--');
y1 = XonOmegaTranspose(U.',V.',omega);   e1 = norm(y1-y0);
y2 = XonOmegaTranspose(U.',V.',Y);       e2 = norm(y2-y0);
y3 = XonOmegaTranspose(U.',V.',I,J);     e3 = norm(y3-y0);
fprintf('\tMethod 1 error: %e\n\tMethod 2 error: %e\n\tMethod 3 error: %e\n',e1,e2,e3);



disp('-- Testing updateSparse --');
A = sprand(100,100,.1) + 1i*sprand(100,100,.1);
[I,J] = find(A); omega = sub2ind( size(A), I, J );
B = zeros(size(A)); B(omega) = 1;  B = complex(B);
B = sparse(B);
updateSparse(B,A(omega));
fprintf('\tDiscrepancy is %f \n', norm(full(B(omega))-A(omega)));

disp('-- Testing smvp --');
% A = sprand(100,100,.1);  % allow A to be complex
x = randn(100,1);
b = smvp(A,x);
fprintf('\tDiscrepancy is %f \n', norm(b - A*x) );
disp('Finished test');