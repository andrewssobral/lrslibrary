%% Tensors
% Tensors are extensions of multidimensional arrays with additional
% operations defined on them. Here we explain the basics of creating and
% working with tensors.
%% Creating a tensor from an array
% The |tensor| command converts a (multidimensional) array to a tensor
% object. 
M = ones(4,3,2); %<-- A 4 x 3 x 2 array.
X = tensor(M) %<-- Convert to a tensor object.
%%
% Optionally, you can specify a different shape for the tensor, so long as
% the input array has the right number of elements.
X = tensor(M,[2 3 4]) %<-- M has 24 elements.
%% Creating a one-dimensional tensor
% The tensor class explicitly supports order-one tensors as well as
% trailing singleton dimensions, but the size must be explicit in the
% constructor. By default, a column array produces a 2-way tensor. 
X = tensor(rand(5,1)) %<-- Creates a 2-way tensor.
%% 
% This is fixed by specifying the size explicitly.
X = tensor(rand(5,1),5) %<-- Creates a 1-way tensor.
%% Specifying trailing singleton dimensions in a tensor
% Likewise, trailing singleton dimensions must be explictly specified.
Y = tensor(rand(4,3,1)) %<-- Creates a 2-way tensor.
%%
Y = tensor(rand(4,3,1),[4 3 1]) %<-- Creates a 3-way tensor.
%%
% Unfortunately, the |whos| command does not report the size of 1D 
% objects correctly (last checked for MATLAB 2006a).
whos X Y %<-- Doesn't report the right size for X!
%% The constituent parts of a tensor
X = tenrand([4 3 2]); %<-- Create data.
X.data %<-- The array.
%%
X.size %<-- The size.
%% Creating a tensor from its constituent parts
Y = tensor(X.data,X.size) %<-- Copies X.
%% Creating an empty tensor
% An empty constructor exists, primarily to support loading previously 
% saved data in MAT-files.
X = tensor %<-- Creates an empty tensor.
%% Use tenone to create a tensor of all ones
X = tenones([3 4 2]) %<-- Creates a 3 x 4 x 2 tensor of ones.
%% Use tenzeros to create a tensor of all zeros
X = tenzeros([1 4 2]) %<-- Creates a 1 x 4 x 2 tensor of zeros.
%% Use tenrand to create a random tensor
X = tenrand([5 4 2]) %<-- Creates a random 5 x 4 x 2 tensor.
%% Use squeeze to remove singleton dimensions from a tensor
squeeze(Y) %<-- Removes singleton dimensions.
%% Use double to convert a tensor to a (multidimensional) array
double(Y) %<-- Converts Y to a standard MATLAB array.
%%
Y.data %<-- Same thing.
%% Use ndims and size to get the size of a tensor
ndims(Y) %<-- Number of dimensions (or ways).
%%
size(Y) %<-- Row vector with the sizes of all dimension. 
%%
size(Y,3) %<-- Size of a single dimension.
%% Subscripted reference for a tensor
X = tenrand([3 4 2 1]); %<-- Create a 3 x 4 x 2 x 1 random tensor.
X(1,1,1,1) %<-- Extract a single element.
%%
% It is possible to extract a subtensor that contains a single element.
% Observe that singleton dimensions are *not* dropped unless they are
% specifically specified, e.g., as above.
X(1,1,1,:) %<-- Produces a tensor of order 1 and size 1.
%%
% In general, specified dimensions are dropped from the result. Here we
% specify the second and third dimension.
X(:,1,1,:) %<-- Produces a tensor of size 3 x 1.
%%
% Moreover, the subtensor is automatically renumbered/resized in the same
% way that MATLAB works for arrays except that singleton dimensions are
% handled explicitly. 
X(1:2,[2 4],1,:) %<-- Produces a tensor of size 2 x 2 x 1.
%%
% It's also possible to extract a list of elements by passing in an array
% of subscripts or a column array of linear indices.
subs = [1,1,1,1; 3,4,2,1]; X(subs) %<-- Extract 2 values by subscript.
%%
inds = [1; 24]; X(inds) %<-- Same thing with linear indices.
%%
% The difference between extracting a subtensor and a list of linear
% indices is ambiguous for 1-dimensional tensors. We can specify 'extract'
% as a second argument whenever we are using a list of subscripts. 
X = tenrand(10); %<-- Create a random tensor.
X([1:6]') %<-- Extract a subtensor.
%%
X([1:6]','extract') %<-- Same thing *but* result is a vector.
%% Subscripted assignment for a tensor 
% We can assign a single element, an entire subtensor, or a list of values
% for a tensor.
X = tenrand([3,4,2]); %<-- Create some data.
X(1,1,1) = 0 %<-- Replaces the (1,1,1) element.
%% 
X(1:2,1:2,1) = ones(2,2) %<-- Replaces a 2 x 2 subtensor.
%% 
X([1 1 1;1 1 2]) = [5;7] %<-- Replaces the (1,1,1) and (1,1,2) elements.
%%
X([1;13]) = [5;7] %<-- Same as above using linear indices.
%%
% It is possible to *grow* the tensor automatically by assigning elements
% outside the original range of the tensor.
X(1,1,3) = 1 %<-- Grows the size of the tensor.
%% Using end for the last array index.
X(end,end,end)  %<-- Same as X(3,4,3).
%%
X(1,1,1:end-1)  %<-- Same as X(1,1,1:2).
%%
% It is also possible to use |end| to index past the end of an array.
X(1,1,end+1) = 5 %<-- Same as X(1,1,4).
%% Use find for subscripts of nonzero elements of a tensor
% The |find| function returns a list of nonzero *subscripts* for a tensor.
% Note that differs from the standard version, which returns linear
% indices.  
X = tensor(floor(3*rand(2,2,2))) %<-- Generate some data.
%%
[S,V] = find(X) %<-- Find all the nonzero subscripts and values.
%%
S = find(X >= 2) %<-- Find subscripts of values >= 2.
%%
V = X(S) %<-- Extract the corresponding values from X.
%% Basic operations (plus, minus, and, or, etc.) on a tensor
% The tensor object supports many basic operations, illustrated here.
A = tensor(floor(3*rand(2,3,2)))
B = tensor(floor(3*rand(2,3,2)))
%%
A & B %<-- Calls and.
%% 
A | B %<-- Calls or.
%%
xor(A,B) %<-- Calls xor.
%%
A==B %<-- Calls eq.
%% 
A~=B %<-- Calls neq.
%% 
A>B %<-- Calls gt.
%%
A>=B %<-- Calls ge.
%%
A<B %<-- Calls lt.
%%
A<=B %<-- Calls le.
%%
~A %<-- Calls not.
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
A.^B %<-- Calls power.
%%
A.^2 %<-- Calls power.
%%
A.\B %<-- Calls ldivide.
%%
A./2 %<-- Calls rdivide.
%%
A./B %<-- Calls rdivide (but beware divides by zero!)
%% Using tenfun for elementwise operations on one or more tensors
% The function |tenfun| applies a specified function to a number of
% tensors. This can be used for any function that is not predefined for
% tensors.
tenfun(@(x)(x+1),A) %<-- Increment every element of A by one.
%%
tenfun(@max,A,B) %<-- Max of A and B, elementwise.
%%
C = tensor(floor(5*rand(2,3,2))) %<-- Create another tensor.
tenfun(@median,A,B,C) %<-- Elementwise means for A, B, and C.
%% Use permute to reorder the modes of a tensor
X = tensor(1:24,[3 4 2]) %<-- Create a tensor.
%%
permute(X,[3 2 1]) %<-- Reverse the modes.
%%
% Permuting a 1-dimensional tensor works correctly.
X = tensor(1:4,4); %<-- Create a 1-way tensor.
permute(X,1) %<-- Call permute with *only* one dimension.
%% Displaying a tensor
% The function |disp| can be used to display a tensor and correctly
% displays very small and large elements.
X = tensor(1:24,[3 4 2]); %<-- Create a 3 x 4 x 2 tensor.
X(:,:,1) = X(:,:,1) * 1e15; %<-- Make the first slice very large.
X(:,:,2) = X(:,:,2) * 1e-15; %<-- Make the second slice very small.
disp(X)
