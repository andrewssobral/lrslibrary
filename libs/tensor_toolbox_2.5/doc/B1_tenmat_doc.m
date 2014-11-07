%% Converting a tensor to a matrix and vice versa
% We show how to convert a tensor to a matrix stored with extra information
% so that it can be converted back to a tensor. Converting to a matrix
% requies an ordered mapping of the tensor indices to the rows and the
% columns of the matrix. 
%% Creating a tenmat (tensor as matrix) object
X = tensor(1:24,[3 2 2 2]) %<-- Create a tensor.
%%
A = tenmat(X,[1 2],[3 4]) %<-- Dims [1 2] map to rows, [3 4] to columns.
%%
B = tenmat(X,[2 1],[3 4]) %<-- Order matters!
%%
C = tenmat(X,[1 2],[4 3]) %<-- Order matters!
%% Creating a tenmat by specifying the dimensions mapped to the rows
% If just the row indices are specified, then the columns are arranged in
% increasing order.
A = tenmat(X,1) %<-- Same as A = tenmat(X,1,2:4)
%% Creating a tenmat by specifying the dimensions mapped to the columns
% Likewise, just the columns can be specified if the 3rd argument is a 't'.
% The rows are arranged in increasing order.
A = tenmat(X, [2 3], 't') %<-- Same as A = tenmat(X,[1 4],[2 3]).
%% Vectorize via tenmat
% All the dimensions can be mapped to the rows or the columnns.
A = tenmat(X,1:4,'t') %<-- Map all the dimensions to the columns
%% Alternative ordering for the columns for mode-n matricization
% Mode-n matricization means that only mode n is mapped to the rows. 
% Different column orderings are available.
A = tenmat(X,2) %<-- By default, columns are ordered as [1 3 4].
%% 
A = tenmat(X,2,[3 1 4]) %<-- Explicit specification.
%%
A = tenmat(X,2,'fc') %<-- Forward cyclic, i.e., [3 4 1].
%%
A = tenmat(X,2,'bc') %<-- Backward cyclic, i.e., [1 4 3].
%% Constituent parts of a tenmat
A.data %<-- The matrix itself.
%%
A.tsize %<-- Size of the original tensor.
%%
A.rdims %<-- Dimensions that were mapped to the rows.
%%
A.cdims %<-- Dimensions that were mapped to the columns.
%% Creating a tenmat from its constituent parts
B = tenmat(A.data,A.rdims,A.cdims,A.tsize) %<-- Recreates A
%% Creating an empty tenmat
B = tenmat %<-- Empty tenmat.
%% Use double to convert a tenmat to a MATLAB matrix
double(A) %<-- converts A to a standard matrix
%% Use tensor to convert a tenmat to a tensor
Y = tensor(A)
%% Use size and tsize for the dimensions of a tenmat
size(A) %<-- Matrix size
tsize(A) %<-- Corresponding tensor size
%% Subscripted reference for a tenmat
A(2,1) %<-- returns the (2,1) element of the matrix.
%% Subscripted assignment for a tenmat
A(1:2,1:2) = ones(2) %<-- Replace part of the matrix.
%% Use end for the last index
A(end,end) %<-- Same as X(2,12)
%% Basic operations for tenmat
norm(A) %<-- Norm of the matrix.
%%
A' %<-- Calls ctranspose (also swaps mapped dimensions).
%%
+A %<-- Calls uplus.
%%
-A %<-- Calls uminus.
%%
A+A %<-- Calls plus.
%%
A-A %<-- Calls minus.
%% Multiplying two tenmats
% It is possible to compute the product of two tenmats and have a result
% that can be converted into a tensor.
B = A * A' %<-- Tenmat that is the product of two tenmats.
%%
tensor(B) %<-- Corresponding tensor.
%% Displaying a tenmat
% Shows the original tensor dimensions, the modes mapped to rows, the modes
% mapped to columns, and the matrix.
disp(A) 
