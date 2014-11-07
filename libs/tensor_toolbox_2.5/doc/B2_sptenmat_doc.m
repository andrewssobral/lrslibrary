%% Converting sparse tensors to matrices and vice versa
% We show how to convert a sptensor to a matrix stored in _coordinate_
% format and with extra information so that it can be converted back to a
% sptensor.

%% Creating a sptenmat (sparse tensor as sparse matrix) object
% A sparse tensor can be converted to a sparse matrix. The matrix, however,
% is not stored as a MATLAB sparse matrix because that format is sometimes
% inefficient for converted sparse tensors. Instead, the row and column
% indices are stored explicitly.
%%
% First, we create a sparse tensor to be converted.
X = sptenrand([10 10 10 10],10) %<-- Generate some data.
%%
% All the same options for tenmat are available as for tenmat.
A = sptenmat(X,1) %<-- Mode-1 matricization.
%%
A = sptenmat(X,[2 3]) %<-- More than one mode is mapped to the columns.
%%
A = sptenmat(X,[2 3],'t') %<-- Specify column dimensions (transpose).
%%
A = sptenmat(X,1:4) %<-- All modes mapped to rows, i.e., vectorize.
%%
A = sptenmat(X,2) %<-- By default, columns are ordered as [1 3 4].
%% 
A = sptenmat(X,2,[3 1 4]) %<-- Explicit column ordering.
%%
A = sptenmat(X,2,'fc') %<-- Foward cyclic.
%%
A = sptenmat(X,2,'bc') %<-- Backward cyclic.
%% Constituent parts of a sptenmat
A.subs %<-- Subscripts of the nonzeros.
%%
A.vals %<-- The corresponding nonzero values.
%%
A.tsize %<-- Size of the original tensor.
%%
A.rdims %<-- Dimensions that were mapped to the rows.
%%
A.cdims %<-- Dimensions that were mapped to the columns.
%% Creating a sptenmat from its constituent parts
B = sptenmat(A.subs,A.vals,A.rdims,A.cdims,A.tsize) %<-- Copies A
%%
B = sptenmat(double(A),A.rdims,A.cdims,A.tsize) %<-- More efficient to pass a matrix.
%% Creating a sptenmat with no nonzeros
A = sptenmat([],[],A.rdims,A.cdims,A.tsize) %<-- An empty sptenmat.
%% Creating an emtpy sptenmat
A = sptenmat %<-- A really empty sptenmat.
%% Use double to convert a sptenmat to a MATLAB sparse matrix
X = sptenrand([10 10 10 10],10); %<-- Create a tensor.
A = sptenmat(X,1) %<-- Convert it to a sptenmat
%%
B = double(A) %<-- Convert it to a MATLAB sparse matrix
%%
whos A B %<-- The storage for B (the sparse matrix) is larger than for A.
%%
C = B'; %<-- Transposing the result fixes the problem.
whos C
%% Use full to convert a sptenmat to a tenmat
B = sptenmat(sptenrand([3 3 3], 3), 1) %<-- Create a sptenmat
%%
C = full(B) %<-- Convert to a tenmat
%% Use sptensor to convert a sptenmat to a sptensor
Y = sptensor(A) %<-- Convert a sptenmat to a sptensor
%% Use size and tsize for the dimensions of a sptenmat
size(A) %<-- Matrix size
tsize(A) %<-- Corresponding tensor size
%% Subscripted reference for a sptenmat
% This is not supported beyond getting the constituent parts.
%% Subscripted assignment for a sptenmat
A(1:2,1:2) = ones(2) %<-- Replace part of the matrix.
%% Use end for the last index
% End is not supported.
%% Basic operations for sptenmat
norm(A) %<-- Norm of the matrix.
%%
+A %<-- Calls uplus.
%%
-A %<-- Calls uminus.
%% Use aatx to efficiently compute A * A' * x for a sptenmat
x = ones(10,1); %<-- Create vector
aatx(A,x) %<-- Compute A * A' * x
%%
double(A) * double(A)' * x %<-- Same as above but less efficient
%% Displaying a tenmat
% Shows the original tensor dimensions, the modes mapped to rows, the modes
% mapped to columns, and the matrix.
disp(A) 
