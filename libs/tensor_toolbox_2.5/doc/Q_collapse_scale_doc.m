%% Collapsing and scaling tensors
% The tensor and sptensor classes support that notion of collapsing and
% scaling dimensions.

%% Examples of collapsing a tensor
X = tenrand([4 3 2]) %<-- Generate some data.
%%
Y = collapse(X,[2 3]) %<-- Sum of entries in each mode-1 slice.
%%
Y = collapse(X,-1) %<-- Same as above.
%%
Z = collapse(X,2) %<-- Sum of entries in each row fiber.
%%
collapse(X,1:3) %<-- Sum of all entries.
%% Alternate accumulation functions for tensor
Y = collapse(X,[1 2],@max) %<-- Max entry in each mode-3 slice.
%%
Z = collapse(X,-3,@mean) %<-- Average entry in each mode-3 slice.
%% Examples of collapsing a sptensor
X = sptenrand([4 3 2],6) %<-- Generate some data.
%%
Y = collapse(X,[2 3]) %<-- Sum of entries in each mode-1 slice.
%%
Y = collapse(X,-1) %<-- Same as above.
%%
Z = collapse(X,2) %<-- Sum of entries in each row fiber.
%%
collapse(X,1:3) %<-- Sum of all entries.
%% Alternate accumulation functions for sptensor
Y = collapse(X,[1 2],@min) %<-- Min *nonzero* entry in each mode-3 slice.
%%
Z = collapse(X,-3,@mean) %<-- Average *nonzero* entry in each mode-3 slice.
%% Scaling a tensor in different modes
X = tenones([3,4,5]); %<-- Generate data 
S = 10 * [1:5]'; Y = scale(X,S,3) %<-- Scale in mode-3
%%
S = tensor(10 * [1:5]',5); Y = scale(X,S,3) %<-- First argument is a tensor.
%%
S = tensor(1:12,[3 4]); Y = scale(X,S,[1 2]) %<-- Scale in two modes.
%%
S = tensor(1:12,[3 4]); Y = scale(X,S,-3) %<-- Same as above.
%%
S = tensor(1:60,[3 4 5]); Y = scale(X,S,1:3) %<-- Scale in every mode.
%%
Y = S .* X %<-- Same as above.

%% Scaling a sptensor in different modes
X = ones(sptenrand([3 4 5], 10)) %<-- Generate data.
%%
S = 10 * [1:5]'; Y = scale(X,S,3) %<-- Scale in one mode.
%%
S = tensor(10 * [1:5]',5); Y = scale(X,S,3) %<-- Same as above.
%%
S = tensor(1:12,[3 4]); Y = scale(X,S,[1 2]) %<-- Scale in two modes.
%%
S = tensor(1:12,[3 4]); Y = scale(X,S,-3) %<-- Same as above.
%%
Z = scale(X,Y,1:3) %<-- Scale by a sparse tensor.
%%
X .* Y %<-- Same as above.