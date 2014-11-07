%% Sparse Tensors
% MATLAB has no native ability to store sparse multidimensional arrays,
% only sparse matrices. Moreover, the compressed sparse column storage
% format for MATLAB sparse matrices is not readily adaptable to sparse
% tensors. Instead, the |sptensor| class stores the data in coordinate
% format.
%% Creating a sptensor 
% A sparse tensor can be created by passing in a list of subscripts and
% values. For example, here we pass in three subscripts and a scalar value.
% The resuling sparse tensor has three nonzero entries, and the size is the
% size of the largest subscript in each dimension.
rand('state',0); %<-- Setup for the script
subs = [1,1,1;1,2,1;3,4,2]; %<-- Subscripts of the nonzeros.
vals = [1; 2; 3]; %<-- The values of the nonzeros.
X = sptensor(subs,vals) %<-- Create a sparse tensor with 3 nonzeros.
%%
X = sptensor(subs,vals,[3 5 2]) %<-- Or, specify the size explicitly.
%%
% Values corresponding to repeated subscripts are summed. Also note that we
% can use a scalar as the second argument.
subs = [1 1 1; 1 1 3; 2 2 2; 4 4 4; 1 1 1; 1 1 1]; %<-- (1,1,1) is repeated.
X = sptensor(subs,2) %<-- Equivalent to X = sptensor(subs,2*ones(6,1)).
%% Specifying the accumulation method for the constructor
% By default, values corresponding to repeated elements are summed.
% However, it is possible to specify other actions to be taken.
X = sptensor(subs,2*ones(6,1),[4 4 4],@max) %<-- Maximum element.
%%
myfun = @(x) sum(x) / 3; %<-- Total sum divided by three.
X = sptensor(subs,2*ones(6,1),[4 4 4],myfun) %<-- Custom accumulation function.
%% Creating a one-dimensional sptensor.
X = sptensor([1;3;5],1,10) %<-- Same as X = sptensor([1;3;5],[1;1;1],1,10).
%%
X = sptenrand(50,5) %<-- A random, sparse, order-1 tensor with 5 nonzeros.
%% Creating an all-zero sptensor 
X = sptensor([],[],[10 10 10]) %<-- Creates an all-zero tensor.
%%
X = sptensor([10 10 10]) %<-- Same as above.
%% Constituent parts of a sptensor
X = sptenrand([40 30 20],5); %<-- Create data.
X.subs %<-- Subscripts of nonzeros.
%%
X.vals %<-- Corresponding nonzero values.
%%
X.size %<-- The size.
%% Creating a sparse tensor from its constituent parts
Y = sptensor(X.subs,X.vals,X.size) %<-- Copies X.
%% Creating an empty sptensor
% An empty constructor exists, primarily to support loads of previously 
% saved data.
Y = sptensor %<-- Create an empty sptensor.
%% Use sptenrand to create a random sptensor
X = sptenrand([10 10 10],0.01) %<-- Create a tensor with 1% nonzeroes.
%%
% It is also posible to specify the precise number of nonzeros rather than
% a percentage.
X = sptenrand([10 10 10],10) %<-- Create a tensor with 10 nonzeros.
%% Use squeeze to remove singleton dimensions from a sptensor
Y = sptensor([1 1 1; 2 1 1], 1, [2 1 1]) %<-- Create a sparse tensor.
squeeze(Y) %<-- Remove singleton dimensions.
%% Use full or tensor to convert a sptensor to a (dense) tensor
X = sptensor([1 1 1; 2 2 2], [1; 1]); %<-- Create a sparse tensor.
Y = full(X) %<-- Convert it to a (dense) tensor.
%%
Y = tensor(X) %<-- Same as above.
%% Use sptensor to convert a (dense) tensor to a sptensor
Z = sptensor(Y) %<-- Convert a tensor to a sptensor.
%% Use double to convert a sptensor to a (dense) multidimensional array
Y = double(X) %<-- Creates a MATLAB array.
%% Use find to extract nonzeros from a tensor and then create a sptensor
% The |find| command can be used to extract specific elements and then
% convert those into a sptensor.
X = tensor(rand(5,4,2),[5 4 2]) %<-- Create a tensor.
S = find(X > 0.9) %<-- Extract subscipts of values greater than 0.9.
V = X(S) %<-- Extract the corresponding values.
Y = sptensor(S,V,[5 4 2]) %<-- Create a new tensor.
%% Use ndims and size to get the size of a sptensor
ndims(Y) %<-- Number of dimensions or modes.
%%
size(Y) %<-- Size of Y.
%% 
size(Y,3) %<-- Size of mode 3 of Y.
%% Use nnz to get the number of nonzeros of a sptensor
nnz(Y) %<-- Number of nonzeros in Y.
%% Subscripted reference for a sptensor
X = sptensor([4,4,4;2,2,1;2,3,2],[3;5;1],[4 4 4]) %<-- Create a sptensor.
%% 
X(1,2,1) %<-- Extract the (1,2,1) element, which is zero.
%%
X(4,4,4) %<-- Extract the (4,4,4) element, which is nonzero.
%% 
X(1:2,2:4,:) %<-- Extract a 2 x 3 x 4 subtensor.
%%
X([1 1 1; 2 2 1]) %<-- Extract elements by subscript.
%%
X([1;6]) %<-- Same as above but with linear indices.
%%
% As with a tensor, subscriped reference may be ambiguous for
% one-dimensional tensors. 
X = sptensor([1;3;5],1,7) %<-- Create a sparse tensor.
%%
X(3) %<-- Fully specified, single elements are always returned as scalars.
%%
X([3;6]) %<-- Returns a subtensor.
%%
X([3;6],'extract') %<-- Same as above *but* returns an array.
%% Subscripted assignment for a sptensor
X = sptensor([30 40 20]) %<-- Create an emtpy 30 x 40 x 20 sptensor.
%% 
X(30,40,20) = 7 %<-- Assign a single element.
%%
X([1,1,1;2,2,2]) = [1;1] %<-- Assign a list of elements.
%%
X(11:20,11:20,11:20) = sptenrand([10,10,10],10) %<-- Assign a subtensor.
%%
X(31,41,21) = 4 %<-- Grows the size of the sptensor.
%%
X(111:120,111:120,111:120) = sptenrand([10,10,10],10) %<-- Grow more.
%% Use end as the last index.
X(end-10:end,end-10:end,end-5:end)  %<-- Same as X(108:118,110:120,115:120)
%% Use elemfun to manipulate the nonzeros of a sptensor
% The function |elemfun| is similar to |spfun| for sparse matrices.
X = sptenrand([10,10,10],3) %<-- Create some data.
%%
Z = elemfun(X, @sqrt) %<-- Square root of every nonzero.
%%
Z = elemfun(X, @(x) x+1) %<-- Use a custom function.
%%
Z = elemfun(X, @(x) x~=0) %<-- Set every nonzero to one.
%%
Z = ones(X) %<-- An easier way to change every nonzero to one.
%% Basic operations (plus, minus, times, etc.) on a sptensor
A = sptensor(tensor(floor(5*rand(2,2,2)))) %<-- Create data.
B = sptensor(tensor(floor(5*rand(2,2,2)))) %<-- Create more data.
%%
+A %<-- Calls uplus.
%%
-A %<-- Calls uminus.
%%
A+B %<-- Calls plus.
%%
A-B %<-- Calls minus.
%%
A.*B %<-- Calls times.
%%
5*A %<-- Calls mtimes.
%%
A./2 %<-- Calls rdivide.
%%
% Elementwise divsion by another sptensor is allowed, but 
% if the sparsity pattern of the denominator should be a
% superset of the numerator.
A./(A+B) %<-- Calls rdivide.
%%
A./B %<-- Uh-oh. Getting a divide by zero.
%% Use permute to reorder the modes of a sptensor
A = sptenrand([30 40 20 1], 5) %<-- Create data.
%%
permute(A,[4 3 2 1]) %<-- Reorder the modes.
%%
% Permute works correctly for a 1-dimensional sptensor.
X = sptenrand(40,4) %<-- Create data.
%%
permute(X,1) %<-- Permute.
%% Displaying a tensor
% The function |disp| handles small and large elements appropriately, as
% well as aligning the indices.
X = sptensor([1 1 1]); %<-- Create an empty sptensor. 
X(1,1,1) = rand(1)*1e15; %<-- Insert a very big element.
X(4,3,2) = rand(1)*1e-15; %<-- Insert a very small element.
X(2,2,2) = rand(1); %<-- Insert a 'normal' element.
disp(X)
