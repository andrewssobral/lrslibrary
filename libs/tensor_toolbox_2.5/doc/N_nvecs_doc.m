%% Generating the leading mode-n vectors
% The leading mode-n vectors are those vectors that span the subspace of
% the mode-n fibers. In other words, the left singular vectors of the
% n-mode matricization of X. 
%% Using nvecs to calculate the leading mode-n vectors
% The |nvecs| command efficient computes the leading n-mode vectors.
rand('state',0);
X = sptenrand([4,3,2],6) %<-- A sparse tensor
%%
nvecs(X,1,2) %<-- The 2 leading mode-1 vectors
%% 
nvecs(X,1,3) % <-- The 3 leading mode-1 vectors
%%
nvecs(full(X),1,3) %<-- The same thing for a dense tensor
%%
X = ktensor({rand(3,2),rand(3,2),rand(2,2)}) %<-- A random ktensor
%%
nvecs(X,2,1) %<-- The 1 leading mode-2 vector
%%
nvecs(full(X),2,1) %<-- Same thing for a dense tensor
%%
X = ttensor(tenrand([2,2,2,2]),{rand(3,2),rand(3,2),rand(2,2),rand(2,2)}); %<-- A random ttensor
%%
nvecs(X,4,2) %<-- The 1 leading mode-2 vector
%%
nvecs(full(X),4,2) %<-- Same thing for a dense tensor
%% Using nvecs for the HOSVD
X = tenrand([4 3 2]) %<-- Generate data
%% 
U1 = nvecs(X,1,4); %<-- Mode 1
U2 = nvecs(X,2,3); %<-- Mode 2
U3 = nvecs(X,3,2); %<-- Mode 3
S = ttm(X,{pinv(U1),pinv(U2),pinv(U3)}); %<-- Core
Y = ttensor(S,{U1,U2,U3}) %<-- HOSVD of X
%%
norm(full(Y) - X) %<-- Reproduces the same result.

%% 
U1 = nvecs(X,1,2); %<-- Mode 1
U2 = nvecs(X,2,2); %<-- Mode 2
U3 = nvecs(X,3,2); %<-- Mode 3
S = ttm(X,{pinv(U1),pinv(U2),pinv(U3)}); %<-- Core
Y = ttensor(S,{U1,U2,U3}) %<-- Rank-(2,2,2) HOSVD approximation of X
%%
100*(1-norm(full(Y)-X)/norm(X)) %<-- Percentage explained by approximation