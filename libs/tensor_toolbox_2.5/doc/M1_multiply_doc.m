%% Multiplying tensors

%% Tensor times vector (ttv for tensor)
% Compute a tensor times a vector (or vectors) in one (or more) modes.
rand('state',0);
X = tenrand([5,3,4,2]); %<-- Create a dense tensor.
A = rand(5,1); B = rand(3,1); C = rand(4,1); D = rand(2,1); %<-- Some vectors.
%%
Y = ttv(X, A, 1) %<-- X times A in mode 1.
%%
Y = ttv(X, {A,B,C,D}, 1) %<-- Same as above.
%%
Y = ttv(X, {A,B,C,D}, [1 2 3 4]) %<-- All-mode multiply produces a scalar.
%%
Y = ttv(X, {D,C,B,A}, [4 3 2 1]) %<-- Same as above.
%%
Y = ttv(X, {A,B,C,D}) %<-- Same as above.
%%
Y = ttv(X, {C,D}, [3 4]) %<-- X times C in mode-3 & D in mode-4.
%%
Y = ttv(X, {A,B,C,D}, [3 4]) %<-- Same as above.
%%
Y = ttv(X, {A,B,D}, [1 2 4]) %<-- 3-way multiplication.
%%
Y = ttv(X, {A,B,C,D}, [1 2 4]) %<-- Same as above.
%%
Y = ttv(X, {A,B,D}, -3) %<-- Same as above.
%%
Y = ttv(X, {A,B,C,D}, -3) %<-- Same as above.
%% Sparse tensor times vector (ttv for sptensor)
% This is the same as in the dense case, except that the result may be
% either dense or sparse (or a scalar).
X = sptenrand([5,3,4,2],5); %<-- Create a sparse tensor.
%%
Y = ttv(X, A, 1) %<-- X times A in mode 1. Result is sparse.
%%
Y = ttv(X, {A,B,C,D}, [1 2 3 4]) %<-- All-mode multiply.
%%
Y = ttv(X, {C,D}, [3 4]) %<-- X times C in mode-3 & D in mode-4. 
%%
Y = ttv(X, {A,B,D}, -3) %<-- 3-way multiplication. Result is *dense*!
%% Kruskal tensor times vector (ttv for ktensor)
% The special structure of a ktensor allows an efficient implementation of
% vector multiplication. The result is a ktensor or a scalar.
X = ktensor([10;1],rand(5,2),rand(3,2),rand(4,2),rand(2,2)); %<-- Ktensor.
Y = ttv(X, A, 1) %<-- X times A in mode 1. Result is a ktensor.
%%
norm(full(Y) - ttv(full(X),A,1)) %<-- Result is the same as dense case.
%%
Y = ttv(X, {A,B,C,D}) %<-- All-mode multiply -- scalar result.
%%
Y = ttv(X, {C,D}, [3 4]) %<-- X times C in mode-3 & D in mode-4.
%%
Y = ttv(X, {A,B,D}, [1 2 4]) %<-- 3-way multiplication.
%% Tucker tensor times vector (ttv for ttensor)
% The special structure of a ttensor allows an efficient implementation of
% vector multiplication. The result is a ttensor or a scalar.
X = ttensor(tenrand([2,2,2,2]),rand(5,2),rand(3,2),rand(4,2),rand(2,2));
Y = ttv(X, A, 1) %<-- X times A in mode 1.
%%
norm(full(Y) - ttv(full(X),A, 1)) %<-- Same as dense case.
%%
Y = ttv(X, {A,B,C,D}, [1 2 3 4]) %<-- All-mode multiply -- scalar result.
%%
Y = ttv(X, {C,D}, [3 4]) %<-- X times C in mode-3 & D in mode-4.
%%
Y = ttv(X, {A,B,D}, [1 2 4]) %<-- 3-way multiplication.
%% Tensor times matrix (ttm for tensor)
% Compute a tensor times a matrix (or matrices) in one (or more) modes.
X = tensor(rand(5,3,4,2));
A = rand(4,5); B = rand(4,3); C = rand(3,4); D = rand(3,2);
%%
Y = ttm(X, A, 1);         %<-- X times A in mode-1.
Y = ttm(X, {A,B,C,D}, 1); %<-- Same as above.
Y = ttm(X, A', 1, 't')    %<-- Same as above.
%%
Y = ttm(X, {A,B,C,D}, [1 2 3 4]); %<-- 4-way mutliply.
Y = ttm(X, {D,C,B,A}, [4 3 2 1]); %<-- Same as above.
Y = ttm(X, {A,B,C,D});            %<-- Same as above.
Y = ttm(X, {A',B',C',D'}, 't')    %<-- Same as above.
%%
Y = ttm(X, {C,D}, [3 4]);    %<-- X times C in mode-3 & D in mode-4
Y = ttm(X, {A,B,C,D}, [3 4]) %<-- Same as above.
%%
Y = ttm(X, {A,B,D}, [1 2 4]);   %<-- 3-way multiply.
Y = ttm(X, {A,B,C,D}, [1 2 4]); %<-- Same as above.
Y = ttm(X, {A,B,D}, -3);        %<-- Same as above.
Y = ttm(X, {A,B,C,D}, -3)       %<-- Same as above.
%% Sparse tensor times matrix (ttm for sptensor)
% It is also possible to multiply an sptensor times a matrix or series of
% matrices. The arguments are the same as for the dense case. The result may
% be dense or sparse, depending on its density. 
X = sptenrand([5 3 4 2],10);
Y = ttm(X, A, 1);         %<-- X times A in mode-1.
Y = ttm(X, {A,B,C,D}, 1); %<-- Same as above.
Y = ttm(X, A', 1, 't')    %<-- Same as above
%%
norm(full(Y) - ttm(full(X),A, 1) ) %<-- Same as dense case.
%%
Y = ttm(X, {A,B,C,D}, [1 2 3 4]); %<-- 4-way multiply.
Y = ttm(X, {D,C,B,A}, [4 3 2 1]); %<-- Same as above.
Y = ttm(X, {A,B,C,D});            %<-- Same as above.
Y = ttm(X, {A',B',C',D'}, 't')    %<-- Same as above.
%%
Y = ttm(X, {C,D}, [3 4]);    %<-- X times C in mode-3 & D in mode-4
Y = ttm(X, {A,B,C,D}, [3 4]) %<-- Same as above.
%%
Y = ttm(X, {A,B,D}, [1 2 4]);   %<-- 3-way multiply.
Y = ttm(X, {A,B,C,D}, [1 2 4]); %<-- Same as above.
Y = ttm(X, {A,B,D}, -3);        %<-- Same as above.
Y = ttm(X, {A,B,C,D}, -3)       %<-- Same as above.
%%
% The result may be dense or sparse. 
X = sptenrand([5 3 4],1);
Y = ttm(X, A, 1) %<-- Sparse result.
%%
X = sptenrand([5 3 4],50);
Y = ttm(X, A, 1) %<-- Dense result.
%%
% Sometimes the product may be too large to reside in memory.  For
% example, try the following:
% X = sptenrand([100 100 100 100], 1e4);
% A = rand(1000,100);
% ttm(X,A,1);  %<-- too large for memory
%% Kruskal tensor times matrix (ttm for ktensor)
% The special structure of a ktensor allows an efficient implementation of
% matrix multiplication. The arguments are the same as for the dense case.
X = ktensor({rand(5,1) rand(3,1) rand(4,1) rand(2,1)});
%%
Y = ttm(X, A, 1);         %<-- X times A in mode-1.
Y = ttm(X, {A,B,C,D}, 1); %<-- Same as above.
Y = ttm(X, A', 1, 't')    %<-- Same as above.
%%
Y = ttm(X, {A,B,C,D}, [1 2 3 4]); %<-- 4-way mutliply.
Y = ttm(X, {D,C,B,A}, [4 3 2 1]); %<-- Same as above.
Y = ttm(X, {A,B,C,D});            %<-- Same as above.
Y = ttm(X, {A',B',C',D'}, 't')    %<-- Same as above.
%%
Y = ttm(X, {C,D}, [3 4]);    %<-- X times C in mode-3 & D in mode-4.
Y = ttm(X, {A,B,C,D}, [3 4]) %<-- Same as above.
%%
Y = ttm(X, {A,B,D}, [1 2 4]);   %<-- 3-way multiply.
Y = ttm(X, {A,B,C,D}, [1 2 4]); %<-- Same as above.
Y = ttm(X, {A,B,D}, -3);        %<-- Same as above.
Y = ttm(X, {A,B,C,D}, -3)       %<-- Same as above.
%% Tucker tensor times matrix (ttm for ttensor)
% The special structure of a ttensor allows an efficient implementation of
% matrix multiplication.
X = ttensor(tensor(rand(2,2,2,2)),{rand(5,2) rand(3,2) rand(4,2) rand(2,2)});
%%
Y = ttm(X, A, 1);         %<-- computes X times A in mode-1.
Y = ttm(X, {A,B,C,D}, 1); %<-- Same as above.
Y = ttm(X, A', 1, 't')    %<-- Same as above.
%%
Y = ttm(X, {A,B,C,D}, [1 2 3 4]); %<-- 4-way multiply.
Y = ttm(X, {D,C,B,A}, [4 3 2 1]); %<-- Same as above.
Y = ttm(X, {A,B,C,D});            %<-- Same as above.
Y = ttm(X, {A',B',C',D'}, 't')    %<-- Same as above.
%%
Y = ttm(X, {C,D}, [3 4]);    %<-- X times C in mode-3 & D in mode-4
Y = ttm(X, {A,B,C,D}, [3 4]) %<-- Same as above.
%%
Y = ttm(X, {A,B,D}, [1 2 4]);   %<-- 3-way multiply
Y = ttm(X, {A,B,C,D}, [1 2 4]); %<-- Same as above.
Y = ttm(X, {A,B,D}, -3);        %<-- Same as above.
Y = ttm(X, {A,B,C,D}, -3)       %<-- Same as above.
%% Tensor times tensor (ttt for tensor)
X = tensor(rand(4,2,3)); Y = tensor(rand(3,4,2));
Z = ttt(X,Y); %<-- Outer product of X and Y.
size(Z)
%%
Z = ttt(X,X,1:3) %<-- Inner product of X with itself.
%%
Z = ttt(X,Y,[1 2 3],[2 3 1]) %<-- Inner product of X & Y.
%%
Z = innerprod(X,Y) %<-- Same as above.
%%
Z = ttt(X,Y,[1 3],[2 1]) %<-- Product of X & Y along specified dims.
%% Sparse tensor times sparse tensor (ttt for sptensor)
X = sptenrand([4 2 3],3); Y = sptenrand([3 4 2],3);
Z = ttt(X,Y) %<--Outer product of X and Y.
%%
norm(full(Z)-ttt(full(X),full(Y))) %<-- Same as dense.
%%
Z = ttt(X,X,1:3) %<-- Inner product of X with itself.
%%
X = sptenrand([2 3],1); Y = sptenrand([3 2],1);
Z = ttt(X, Y) %<-- Sparse result.
%%
X = sptenrand([2 3],20); Y = sptenrand([3 2],20);
Z = ttt(X, Y) %<-- Dense result.
%%
Z = ttt(X,Y,[1 2],[2 1]) %<-- inner product of X & Y
%% Inner product (innerprod)
% The function |innerprod| efficiently computes the inner product
% between two tensors X and Y.  The code does this efficiently
% depending on what types of tensors X and Y.
X = tensor(rand(2,2,2))
Y = ktensor({rand(2,2),rand(2,2),rand(2,2)})
%%
z = innerprod(X,Y)
%% Contraction on tensors (contract for tensor)
% The function |contract| sums the entries of X along dimensions I and
% J.  Contraction is a generalization of matrix trace. In other words,
% the trace is performed along the two-dimensional slices defined by
% dimensions I and J. It is possible to implement tensor
% multiplication as an outer product followed by a contraction.
X = sptenrand([4 3 2],5); 
Y = sptenrand([3 2 4],5);
%%
Z1 = ttt(X,Y,1,3); %<-- Normal tensor multiplication
%%
Z2 = contract(ttt(X,Y),1,6); %<-- Outer product + contract
%%
norm(Z1-Z2) %<-- Should be zero
%%
% Using |contract| on either sparse or dense tensors gives the same
% result
X = sptenrand([4 2 3 4],20); 
Z1 = contract(X,1,4)        % sparse version of contract
%%
Z2 = contract(full(X),1,4)  % dense version of contract
%%
norm(full(Z1) - Z2) %<-- Should be zero
%%
% The result may be dense or sparse, depending on its density. 
X = sptenrand([4 2 3 4],8); 
Y = contract(X,1,4) %<-- should be sparse
%%
X = sptenrand([4 2 3 4],80);
Y = contract(X,1,4) %<-- should be dense
%% Relationships among ttv, ttm, and ttt 
% The three "tensor times ___" functions (|ttv|, |ttm|, |ttt|) all perform
% specialized calculations, but they are all related to some degree.
% Here are several relationships among them:
%%
X = tensor(rand(4,3,2)); 
A = rand(4,1);
%%
% Tensor times vector gives a 3 x 2 result
Y1 = ttv(X,A,1)
%%
% When |ttm| is used with the transpose option, the result is almost
% the same as |ttv|
Y2 = ttm(X,A,1,'t')
%%
% We can use |squeeze| to remove the singleton dimension left over
% from |ttm| to give the same answer as |ttv|
squeeze(Y2)
%%
% Tensor outer product may be used in conjuction with contract to
% produce the result of |ttm|.  Please note that this is more expensive
% than using |ttm|.
Z = ttt(tensor(A),X);
size(Z)
%%
Y3 = contract(Z,1,3)
%%
% Finally, use |squeeze| to remove the singleton dimension to get
% the same result as |ttv|.
squeeze(Y3)
%% Frobenius norm of a tensor
% The Frobenius norm of any type of tensor may be computed with the 
% function |norm|.  Each class is optimized to calculate the norm
% in the most efficient manner.
X = sptenrand([4 3 2],5)
norm(X)
norm(full(X))
%%
X = ktensor({rand(4,2),rand(3,2)})
norm(X)
%%
X = ttensor(tensor(rand(2,2)),{rand(4,2),rand(3,2)})
norm(X)
