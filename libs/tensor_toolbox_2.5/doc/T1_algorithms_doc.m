%% ALS optimization for CP and Tucker tensor decompositions

%% Alternating least squares for PARAFAC/CANDECOMP
% The function |parafac_als| computes an estimate of the best rank-R
% PARAFAC model of a tensor X using an alternating least-squares
% algorithm.  The input X can be a tensor, sptensor, ktensor, or
% ttensor. The result P is a ktensor.
rand('state',0);
X = sptenrand([5 4 3], 10)
%%
P = parafac_als(X,2)
%%
P = parafac_als(X,2,struct('dimorder',[3 2 1]))
%%
P = parafac_als(X,2,struct('dimorder',[3 2 1],'init','nvecs'))
%%
U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of P
P = parafac_als(X,2,struct('dimorder',[3 2 1],'init',{U0}))
%% Alternating least squares for Tucker model 
% The function |tucker_als| computes the best rank(R1,R2,..,Rn)
% approximation of tensor X, according to the specified dimensions in
% vector R.  The input X can be a tensor, sptensor, ktensor, or
% ttensor.  The result returned in T is a ttensor.
X = sptenrand([5 4 3], 10)
%%
T = tucker_als(X,2)        %<-- best rank(2,2,2) approximation 
%%
T = tucker_als(X,[2 2 1])  %<-- best rank(2,2,1) approximation 
%%
T = tucker_als(X,2,struct('dimorder',[3 2 1]))
%%
T = tucker_als(X,2,struct('dimorder',[3 2 1],'init','eigs'))
%%
U0 = {rand(5,2),rand(4,2),[]}; %<-- Initial guess for factors of T
T = tucker_als(X,2,struct('dimorder',[3 2 1],'init',{U0}))
