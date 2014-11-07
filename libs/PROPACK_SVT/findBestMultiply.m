function [forwardType,transposeType,details] = findBestMultiply(A,T)
% [forwardType,transposeType] = findBestMultiply(A)
%   does a speed test to determine what the fastest routine
%   for sparse matrix multiplication is, for the matrix A.
% ... = findBestMultiply(A,T)
%   tries to make each test last T seconds (there are a 6 tests total)
%
% The fastest routine will depend on the verson of Matlab, the
% computer OS, the cpu type (dual core? etc.), and the particular
% size and structure of the matrix.
%
% Stephen Becker, srbecker@caltech.edu, 3/14/09

% 5/13/09, smvp now handles complex data, so modify this accordingly

if nargin < 2, T = 1; end

At = A';
x = randn(size(A,2),1);
y = randn(size(A,1),1);
if ~isreal(A)
    x = x + 1i*randn(size(A,2),1);
    y = y + 1i*randn(size(A,1),1);
end
    

% how many times to repeat the test?
tic; A'*y; t = toc;

% if we want each test to last say, T, seconds, then...
nIter = max( 1, round(T/t) );


% -- for the forward multiply --

% -- Method 1: use standard MATLAB
tic; for i = 1:nIter, A*x; end; t1 = toc;

% -- Method 2: use transpose trick (works best in new versions)
tic; for i = 1:nIter, At'*x; end; t2 = toc;

% -- Method 3: simple mex file
try
    tic; for i = 1:nIter, smvp(A,x); end; t3 = toc;
catch
    l = lasterror;
    fprintf('Error thrown: %s\n',l.message);
    disp('findBestMultiply: Sorry, this probably means the smvp mex file\n\t is not installed');
    disp('Try compiling it with:  mex -O smvp.c  in the "private" subdirectory');
    t3 = Inf;
end

[m,forwardType] = min( [t1,t2,t3] );
if nargout > 2, details = [t1,t2,t3]; end


% -- for the transpose multiply --

% -- Method 1: use standard MATLAB
tic; for i = 1:nIter, At*y; end; t1 = toc;

% -- Method 2: use transpose trick (works best in new versions)
tic; for i = 1:nIter, A'*y; end; t2 = toc;

% -- Method 3: simple mex file
try
    tic; for i = 1:nIter, smvp(At,y); end; t3 = toc;
catch
%     disp('Sorry, smvp mex file not installed');
    t3 = Inf;
end

[m,transposeType] = min( [t1,t2,t3] );
if nargout > 2, details = [details,[t1,t2,t3] ]; end