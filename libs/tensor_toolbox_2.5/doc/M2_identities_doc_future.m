%% Identities and relationships of tensors
% There are many mathematical relationships, identities, and
% connections among tensors.  These identities are presented here and
% show the versatility of the Tensor Toolbox.
% The propositions indicated below are references to the following
% report:
%
% T.G. Kolda, "Multilinear operators for higher-order decompositions",
% Tech. Rep. SAND2006-2081, Sandia National Laboratories, April 2006,
% http://csmr.ca.sandia.gov/~tgkolda/pubs/index.html#SAND2006-2081.

%% N-mode product properties
% Create some data.
Y = tenrand([4 3 2]);
A = rand(3,4);
B = rand(3,3);
%%
% Prop 3.4(a): The order of the multiplication in different modes is irrelevant. 
%
% $$(Y \times_1 A) \times_2 B = (Y \times_2 B) \times_1 A$$
%
X1 = ttm( ttm(Y,A,1), B, 2); %<-- Y x_1 A x_2 B
X2 = ttm( ttm(Y,B,2), A, 1); %<-- Y x_2 B x_1 A
norm(X1 - X2) %<-- difference is zero
%% N-mode product and matricization
% Generate some data to work with.
Y = tenrand([5 4 3]);
A = rand(4,5); B = rand(3,4); C = rand(2,3); U = {A,B,C};
%%
% Prop. 3.7a: N-mode multiplication can be expressed in terms of matricized
% tensors.
%
% $$X = Y \times_n U \Leftrightarrow  X_{(n)} = UY_{(n)} $$
% 
for n = 1:ndims(Y)
  X = ttm(Y,U,n); %<-- X = Y x_n U{n}
  Xn = U{n} * tenmat(Y,n); %<-- Xn = U{n} * Yn
  norm(tenmat(X,n) - Xn)  % <-- should be zero
end
%%
% Prop. 3.7b: We can do matricizations in various ways and still be
% equivalent.
X = ttm(Y,U); %<-- X = Y x_1 A x_2 B x_3 C
Xm1 = kron(B,A)*tenmat(Y,[1 2])*C';  %<-- Kronecker product version
Xm2 = tenmat(X,[1 2]); %<-- Matriczed version
norm(Xm1 - Xm2)  % <-- should be zero
Xm1 = B * tenmat(Y,2,[3 1]) * kron(A,C)'; %<-- Kronecker product version
Xm2 = tenmat(X,2,[3 1]); %<-- Matricized version
norm(Xm1 - Xm2) % <-- should be zero
Xm1 = tenmat(Y,[],[1 2 3]) * kron(kron(C,B),A)'; %<-- Vectorized via Kronecker
Xm2 = tenmat(X,[],[1 2 3]); %<-- Vectorized via matricize
norm(Xm1 - Xm2)

%% Norm of difference between two tensors
% Prop. 3.9: For tensors X and Y, we have:
%
% $$\|X-Y\|^2 = \|X\|^2 + \|Y\|^2 - 2<X,Y> $$
%
X = tenrand([5 4 3]); Y = tenrand([5 4 3]);
% The following 2 results should be equal
norm(X-Y)
sqrt(norm(X)^2 - 2*innerprod(X,Y) + norm(Y)^2)
%% 
% This relationship makes it more convenient to compare the norm of
% the difference between two different tensor objects.  Imagine if we
% have a |sptensor| and a |ktensor| and we want the norm of the
% difference, which may be needed to check for convergence, for
% example, but which is very expensive to convert to a full (dense)
% tensor.  Because |innerprod| and |norm| are defined for all types of
% tensor objects, this is a handy formula.
X = sptensor(X);
Y = ktensor({[1:5]',[1:4]',[1:3]'});
% The following 2 results should be equal
norm(full(X)-full(Y))
sqrt(norm(X)^2 - 2*innerprod(X,Y) + norm(Y)^2)

%% Tucker tensor properties
% The properties of the Tucker operator follow directly from the
% properties of n-mode multiplication.

% Initialize data
Y = tensor(1:24,[4 3 2]);
A1 = reshape(1:20,[5 4]);
A2 = reshape(1:12,[4 3]);
A3 = reshape(1:6,[3 2]);
A = {A1,A2,A3};
B1 = reshape(1:20,[4 5]);
B2 = reshape(1:12,[3 4]);
B3 = reshape(1:6,[2 3]);
B = {B1,B2,B3};
%%
% Proposition 4.2a
X = ttensor(ttensor(Y,A),B)
%%
AB = {B1*A1, B2*A2, B3*A3};
Y = ttensor(Y,AB)
%%
norm(full(X)-full(Y))  %<-- should be zero
%%
% Proposition 4.2b
Y = tensor(1:24,[4 3 2]);
X = ttensor(Y,A);
Apinv = {pinv(A1),pinv(A2),pinv(A3)};
Y2 = ttensor(full(X),Apinv);
norm(full(Y)-full(Y2))  %<-- should be zero
%%
% Proposition 4.2c
Y = tensor(1:24,[4 3 2]);
rand('state',0);
Q1 = orth(rand(5,4));
Q2 = orth(rand(4,3));
Q3 = orth(rand(3,2));
Q = {Q1,Q2,Q3};
X = ttensor(Y,Q)
%%
Qt = {Q1',Q2',Q3'};
Y2 = ttensor(full(X),Qt)
norm(full(Y)-full(Y2))  %<-- should be zero

%% Tucker operator and matricized tensors
% The Tucker operator also has various epressions in terms of
% matricized tensors and the Kronecker product.
% Proposition 4.3a
Y = tensor(1:24,[4 3 2]);
A1 = reshape(1:20,[5 4]);
A2 = reshape(1:12,[4 3]);
A3 = reshape(1:6,[3 2]);
A = {A1,A2,A3};
X = ttensor(Y,A)
for n = 1:ndims(Y)
  rdims = n;
  cdims = setdiff(1:ndims(Y),rdims);
  Xn = A{n} * tenmat(Y,rdims,cdims) * kron(A{cdims(2)}, A{cdims(1)})';
  norm(tenmat(full(X),rdims,cdims) - Xn)  % <-- should be zero
end

%% Orthogonalization of Tucker factors
% Proposition 4.4
Y = tensor(1:24,[4 3 2]);
A1 = rand(5,4);
A2 = rand(4,3);
A3 = rand(3,2);
A = {A1,A2,A3};
X = ttensor(Y,A)
%%
[Q1,R1] = qr(A1);
[Q2,R2] = qr(A2);
[Q3,R3] = qr(A3);
R = {R1,R2,R3};
Z = ttensor(Y,R);
norm(X) - norm(Z)  %<-- should be zero

%% Kruskal operator properties
% Proposition 5.2
A1 = reshape(1:10,[5 2]);
A2 = reshape(1:8,[4 2]);
A3 = reshape(1:6,[3 2]);
K = ktensor({A1,A2,A3});
B1 = reshape(1:20,[4 5]);
B2 = reshape(1:12,[3 4]);
B3 = reshape(1:6,[2 3]);
X = ttensor(K,{B1,B2,B3})

Y = ktensor({B1*A1, B2*A2, B3*A3});
norm(full(X) - full(Y))  %<-- should be zero

%%
% Proposition 5.3a (second part)
A1 = reshape(1:10,[5 2]);
A2 = reshape(1:8,[4 2]);
A3 = reshape(1:6,[3 2]);
A = {A1,A2,A3};
X = ktensor(A);
rdims = 1:ndims(X);
Z = double(tenmat(full(X), rdims, []));
Xn = khatrirao(A{rdims},'r') * ones(length(X.lambda),1);
norm(Z - Xn)  % <-- should be zero
%%
cdims = 1:ndims(X);
Z = double(tenmat(full(X), [], cdims));
Xn = ones(length(X.lambda),1)' * khatrirao(A{cdims},'r')';
norm(Z - Xn)  % <-- should be zero
%%
% Proposition 5.3b
A1 = reshape(1:10,[5 2]);
A2 = reshape(1:8,[4 2]);
A3 = reshape(1:6,[3 2]);
A = {A1,A2,A3};
X = ktensor(A);
for n = 1:ndims(X)
  rdims = n;
  cdims = setdiff(1:ndims(X),rdims);
  Xn = khatrirao(A{rdims}) * khatrirao(A{cdims},'r')';
  Z = double(tenmat(full(X),rdims,cdims));
  norm(Z - Xn)  % <-- should be zero
end
%%
% Proposition 5.3a (first part)
X = ktensor(A);
for n = 1:ndims(X)
  cdims = n;
  rdims = setdiff(1:ndims(X),cdims);
  Xn = khatrirao(A{rdims},'r') * khatrirao(A{cdims})';
  Z = double(tenmat(full(X),rdims,cdims));
  norm(Z - Xn)  % <-- should be zero
end

%% Norm of Kruskal operator
% The norm of a ktensor has a special form because it can be
% reduced to summing the entries of the Hadamard product of N
% matrices of size R x R.
% Proposition 5.4
A1 = reshape(1:10,[5 2]);
A2 = reshape(1:8,[4 2]);
A3 = reshape(1:6,[3 2]);
A = {A1,A2,A3};
X = ktensor(A);
M = ones(size(A{1},2), size(A{1},2));
for i = 1:numel(A)
  M = M .* (A{i}'*A{i});
end
norm(X) - sqrt(sum(M(:)))  %<-- should be zero

%% Inner product of Kruskal operator with a tensor
% The inner product of a ktensor with a tensor yields
% Proposition 5.5
X = tensor(1:60,[5 4 3]);
A1 = reshape(1:10,[5 2]);
A2 = reshape(2:9,[4 2]);
A3 = reshape(3:8,[3 2]);
A = {A1,A2,A3};
K = ktensor(A);
v = khatrirao(A,'r') * ones(size(A{1},2),1);
% The following 2 results should be equal
double(tenmat(X,1:ndims(X),[]))' * v 
innerprod(X,K)
