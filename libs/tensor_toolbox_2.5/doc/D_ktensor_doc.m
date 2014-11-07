%% Kruskal tensors
% Kruskal format is a decomposition of a tensor X as the sum of the outer
% products a the columns of matrices. For example, we might write
% 
% $${\mathcal X} = \sum_r a_r \circ b_r \circ c_r$$
% 
% where a subscript denotes column index and a circle denotes outer
% product. In other words, the tensor X is built frm the columns of the
% matrices A,B, and C. It's often helpful to explicitly specify a weight
% for each outer product, which we do here:
% 
% $${\mathcal X} = \sum_r \lambda_r \; a_r \circ b_r \circ c_r$$
% 
% The |ktensor| class stores the components of the tensor X and can perform
% many operations, e.g., |ttm|, without explicitly forming the tensor X. 

%% Kruskal tensor format via ktensor
% Kruskal format stores a tensor as a sum of rank-1 outer products. For
% example, consider a tensor of the following form.
%  
% $$X = a_1 \circ b_1 \circ c_1 + a_2 \circ b_2 \circ c_2$$
% 
% This can be stored in Kruskal form as follows.
rand('state',0);
A = rand(4,2); %<-- First column is a_1, second is a_2.
B = rand(3,2); %<-- Likewise for B.
C = rand(2,2); %<-- Likewise for C.
X = ktensor({A,B,C}) %<-- Create the ktensor.
%%
% For Kruskal format, there can be any number of matrices, but every matrix
% must have the same number of columns. The number of rows can vary.
Y = ktensor({rand(4,1),rand(2,1),rand(3,1)}) %<-- Another ktensor.
%% Specifying weights in a ktensor
% Weights for each rank-1 tensor can be specified by passing in a 
% column vector.  For example, 
%  
% $$X = \lambda_1 \; a_1 \circ b_1 \circ c_1 + \lambda_2 \; a_2 \circ b_2 \circ c_2$$
% 
lambda = [5.0; 0.25]; %<-- Weights for each factor.
X = ktensor(lambda,{A,B,C}) %<-- Create the ktensor.
%% Creating a one-dimensional ktensor
Y = ktensor({rand(4,5)}) %<-- A one-dimensional ktensor.
%% Constituent parts of a ktensor
X.lambda %<-- Weights or multipliers.
%%
X.U %<-- Cell array of matrices.
%% Creating a ktensor from its constituent parts
Y = ktensor(X.lambda,X.U) %<-- Recreate X.
%% Creating an empty ktensor
Z = ktensor %<-- Empty ktensor.
%% Use full or tensor to convert a ktensor to a tensor
full(X) %<-- Converts to a tensor. 
%% 
tensor(X) %<-- Same as above.
%% Use double to convert a ktensor to a multidimensional array
double(X) %<-- Converts to an array.
%% Use tendiag or sptendiag to convert a ktensor to a ttensor.
% A ktensor can be regarded as a ttensor with a diagonal core.
R = length(X.lambda);  %<-- Number of factors in X.
core = tendiag(X.lambda, repmat(R,1,ndims(X))); %<-- Create a diagonal core.
Y = ttensor(core, X.u) %<-- Assemble the ttensor.
%%
norm(full(X)-full(Y)) %<-- They are the same.
%% 
core = sptendiag(X.lambda, repmat(R,1,ndims(X))); %<-- Sparse diagonal core.
Y = ttensor(core, X.u) %<-- Assemble the ttensor
%%
norm(full(X)-full(Y)) %<-- They are the same.
%% Use ndims and size for the dimensions of a ktensor
ndims(X) %<-- Number of dimensions.
%%
size(X) %<-- Row vector of the sizes.
%%
size(X,2) %<-- Size of the 2nd mode.
%% Subscripted reference for a ktensor
X(1,1,1) %<-- Assemble the (1,1,1) element (requires computation).
%%
X.lambda(2) %<-- Weight of 2nd factor.
%%
X.U{2} %<-- Extract a matrix.
%%
X{2} %<-- Same as above.
%% Subscripted assignment for a ktensor
X.lambda = ones(size(X.lambda)) %<-- Insert new multipliers.
%%
X.lambda(1) = 7 %<-- Change a single element of lambda.
%%
X{3}(1:2,1) = [1;1] %<-- Change the matrix for mode 3.
%% Use end for the last array index.
X(3:end,1,1)  %<-- Calculated X(3,1,1) and X((4,1,1).
%%
X(1,1,1:end-1)  %<-- Calculates X(1,1,1).
%%
X{end}  %<-- Or use inside of curly braces. This is X{3}.
%% Adding and subtracting ktensors
% Adding two ktensors is the same as concatenating the matrices
X = ktensor({rand(4,2),rand(2,2),rand(3,2)}) %<-- Data.
Y = ktensor({rand(4,2),rand(2,2),rand(3,2)}) %<-- More data.
%%
Z = X + Y %<-- Concatenates the factor matrices.
%%
Z = X - Y %<-- Concatenates as with plus, but changes the weights.
%%
norm( full(Z) - (full(X)-full(Y)) ) %<-- Should be zero.
%% Basic operations with a ktensor
+X %<-- Calls uplus.
%%
-X %<-- Calls uminus.
%%
5*X %<-- Calls mtimes.
%% Use permute to reorder the modes of a ktensor
permute(X,[2 3 1]) %<-- Reorders modes of X
%% Use arrange to normalize the factors of a ktensor
% The function |arrange| normalizes the columns of the factors and then
% arranges the rank-one pieces in decreasing order of size.
X = ktensor({rand(3,2),rand(4,2),rand(2,2)})  % <-- Unit weights.
%%
arrange(X) %<-- Normalized a rearranged.
%% Use fixsigns for sign indeterminacies in a ktensor
% The largest magnitude entry for each factor is changed to be
% positive provided that we can flip the signs of _pairs_ of vectors in
% that rank-1 component.
Y = X;
Y.u{1}(:,1) = -Y.u{1}(:,1);  % switch the sign on a pair of columns
Y.u{2}(:,1) = -Y.u{2}(:,1)
%%
fixsigns(Y)
%% Use ktensor to store the 'skinny' SVD of a matrix
A = rand(4,3) %<-- A random matrix.
%%
[U,S,V] = svd(A,0); %<-- Compute the SVD.
X = ktensor(diag(S),{U,V}) %<-- Store the SVD as a ktensor.
%%
double(X) %<-- Reassemble the original matrix.
%% Displaying a ktensor
disp(X) %<-- Displays the vector lambda and each factor matrix.
%% Displaying data
% The |datadisp| function allows the user to associate meaning to the modes
% and display those modes with the most meaning (i.e., corresponding to the
% largest values). 
X = ktensor({[0.8 0.1 1e-10]',[1e-5 2 3 1e-4]',[0.5 0.5]'}); %<-- Create tensor.
X = arrange(X) %<-- Normalize the factors.
%%
labelsDim1 = {'one','two','three'}; %<-- Labels for mode 1.
labelsDim2 = {'A','B','C','D'}; %<-- Labels for mode 2.
labelsDim3 = {'on','off'}; %<-- Labels for mode 3.
datadisp(X,{labelsDim1,labelsDim2,labelsDim3}) %<-- Display.