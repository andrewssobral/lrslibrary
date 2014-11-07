%% Tucker Tensors
% Tucker format is a decomposition of a tensor X as the product of a core
% tensor G and matrices (e.g., A,B,C) in each dimension. In other words, a
% tensor X is expressed as:
% 
% $${\mathcal X} = {\mathcal G} \times_1 A \times_2 B \times_2 C$$
% 
% In MATLAB notation, |X=ttm(G,{A,B,C})|. The |ttensor| class stores the
% components of the tensor X and can perform many operations, e.g., |ttm|,
% without explicitly forming the tensor X.
%% Creating a ttensor with a tensor core
core = tensor(rand(3,2,1),[3 2 1]); %<-- The core tensor.
U = {rand(5,3), rand(4,2), rand(3,1)}; %<-- The matrices.
X = ttensor(core,U) %<-- Create the ttensor.
%% Alternate core formats: sptensor, ktensor, or ttensor
core1 = sptenrand([3 2 1],3); %<-- Create a 3 x 2 x 1 sptensor.
Y = ttensor(core1,U) %<-- Core is a sptensor.
%%
V = {rand(3,2),rand(2,2),rand(1,2)}; %<-- Create some random matrices.
core2 = ktensor(V); %<-- Create a 3 x 2 x 1 ktensor.
Y = ttensor(core2,U) %<-- Core is a ktensor.
%% 
core3 = ttensor(tensor(1:8,[2 2 2]),V); %<-- Create a 3 x 2 x 1 ttensor.
Y = ttensor(core3,U) %<-- Core is a ttensor.
%% Creating a one-dimensional ttensor
Z = ttensor(tensor(rand(2,1),2), rand(4,2)) %<-- One-dimensional ttensor.
%% Constituent parts of a ttensor
X.core %<-- Core tensor.
%%
X.U %<-- Cell array of matrices.
%% Creating a ttensor from its constituent parts
Y = ttensor(X.core,X.U) %<-- Recreate a tensor from its parts.
%% Creating an empty ttensor.
X = ttensor %<-- empty ttensor
%% Use full or tensor to convert a ttensor to a tensor
X = ttensor(core,U) %<-- Create a tensor
%%
full(X) %<-- Converts to a tensor.
%%
tensor(X) %<-- Also converts to a tensor.
%% Use double to convert a ttensor to a (multidimensional) array
double(X) %<-- Converts to a MATLAB array
%% Use ndims and size to get the size of a ttensor
ndims(X) %<-- Number of dimensions.
%%
size(X) %<-- Row vector of the sizes.
%%
size(X,2) %<-- Size of the 2nd mode.
%% Subscripted reference to a ttensor
X.core(1,1,1) %<-- Access an element of the core.
%%
X.U{2} %<-- Extract a matrix.
%%
X{2} %<-- Same as above.
%% Subscripted assignment for a ttensor
X.core = tenones(size(X.core)) %<-- Insert a new core.
%%
X.core(2,2,1) = 7 %<-- Change a single element.
%%
X{3}(1:2,1) = [1;1] %<-- Change the matrix for mode 3.
%% Using end for last index
X{end}  %<-- The same as X{3}.
%% Basic operations (uplus, uminus, mtimes) for a ttensor.
X = ttensor(tenrand([2 2 2]),{rand(3,2),rand(1,2),rand(2,2)}) %<-- Data.
+X %<-- Calls uplus.
%%
-X %<-- Calls uminus.
%%
5*X %<-- Calls mtimes.
%% Use permute to reorder the modes of a ttensor
permute(X,[3 2 1]) %<-- Reverses the modes of X
%% Displaying a ttensor
% The tensor displays by displaying the core and each of the component
% matrices.
disp(X) %<-- Prints out the ttensor.


