function tucker_me_test
%TUCKER_ME_TEST Very simple tests of tucker_me.
%   Code by Tamara Kolda and Jimeng Sun, 2008. 
%
%   Based on the paper:
%   T. G. Kolda and J. Sun. Scalable Tensor Decompositions for Multi-aspect
%   Data Mining. In: ICDM 2008: Proceedings of the 8th IEEE International
%   Conference on Data Mining, December 2008.  

% $Id: tucker_me_test.m,v 1.2 2009/07/08 00:03:52 tgkolda Exp $


%% Set up
csz = [3 6 9 12];
tsz = [50 50 50 50];
X = tucker_me_test_gendata(csz, tsz, .001);

%%
fprintf('-----------------------------------------------------------\n');
fprintf('%-20s | %-7s | %-10s |%-10s \n','Method','Time(s)','Error', 'Memory ratio');
fprintf('-----------------------------------------------------------\n');
%%
tic
[T,Uinit] = tucker_als(X,csz);
t = toc;
fprintf('%-20s | %7.2f | %10.3e \n', 'tucker_als (standard)', t, 1-norm(T)/norm(X));
 
% Make sure that they all use the same initial guess
opts.init = Uinit;

tic
[T, mem_orig] = tucker_me(X,csz, 0, opts);
t = toc;
fprintf('%-20s | %7.2f | %10.3e |%f \n', 'tucker_me (standard)', t, 1-norm(T)/norm(X), 1);
 

%%
tic
[T,mem] = tucker_me(X,csz,1,opts);
t = toc;
fprintf('%-20s | %7.2f | %10.3e|%f  \n', 'tucker_me (slicewise)', t, 1-norm(T)/norm(X), mem/mem_orig);

%%
tic
[T,mem] = tucker_me(X,csz,2,opts);
t = toc;
fprintf('%-20s | %7.2f | %10.3e |%f \n', 'tucker_me (fiberwise)', t, 1-norm(T)/norm(X), mem/mem_orig);

%% Really slow! 
tic
[T,mem] =  tucker_me(X,csz,3,opts);
t = toc;
fprintf('%-20s | %7.2f | %10.3e |%f \n', 'tucker_me (elementwise)', t, 1-norm(T)/norm(X), mem/mem_orig);

%----------------------------------------------------------------------
function [X,G,U] = tucker_me_test_gendata(csz,tsz,pnz)
%TUCKER_ME_TEST_GENDATA Generate sparse array for tucker_me tests.
%
%   X = GENDATA(CSZ,TSZ,PNZ) generates a tensor as follows. 1) Randomly
%   generates a core tensor, G, of size CSZ. 2) Expand it into a tensor of
%   size  TSZ by multiplying it by appropriately sized random matrices in
%   each mode. 3) Save only the largest nonzeros, as specified by PNZ, the
%   percentage of nonzeros to be saved.
%
%   [X,G,U] = GENDATE(...) also returns the core tensor and matrices that
%   were use to generate X.
%
%   Example
%   X = gendata([2 2 2], [10 10 10], .10) - generates a tensor of size 10 x
%   10 x 10 that came from a core of size 2 x 2 x 2 but where the smallest
%   entries were deleted.
%
%   See also TUCKER_ME_TEST.
%
%   Code by Tamara Kolda and Jimeng Sun, 2008. 
%
%   Based on the paper:
%   T. G. Kolda and J. Sun. Scalable Tensor Decompositions for Multi-aspect
%   Data Mining. In: ICDM 2008: Proceedings of the 8th IEEE International
%   Conference on Data Mining, December 2008.  


% Number of tensor dimensions
N = length(csz);

% Generate random core tensor
G = tenrand(csz);

% Generate random U-matrices
for n = 1:N
    U{n} = rand( [tsz(n), csz(n)] );
end

% Create full tensor
Xfull = ttm(G,U);

% Sparsify
[v,idx] = sort(Xfull(:),'descend');
stop = round( pnz * length(idx) );
v(stop:end) = 0;
X = sptensor( reshape(v,size(Xfull)) );

 


