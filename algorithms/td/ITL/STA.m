function [Tnew, S] = STA(Xnew, R, T, S, alpha, samplingPercent)
%Streaming tensor decomposition STA 
%   
%   [Tnew, S] = STA(Xnew, R, T, alpha, samplingPercent) approximately updates tensor PCA,
%   according to new tensor Xnew, old tucker decomposition T and he
%   specified dimensions in vector R. The input Xnew is a tensor or sptensor.
%   the result returned in T is a new tucker tensor for Xnew, S has the
%   energy vector along each mode.
%
%   Examples:
%         % construct two low rank tensor A and B
%         >>A = full(ttensor(tenrand([2,3,4]),{rand(10,2),rand(30,3), rand(40,4)})); 
%         >>B = full(ttensor(tenrand([2,3,4]),{rand(10,2),rand(30,3), rand(40,4)})); 
%         % create a tensor sequence/stream D with 5 copies of A and B
%         % respectively
%         >>D = {A,A,A,A,A,B,B,B,B,B};
%         % compute STA with no forgetting factor
%         >>[T,S] = STA(D{1},[2,3,4]);
%         >>for i = 2:10
%               [T,S] = STA(D{i},[2,3,4],T,S);
%               err = norm(full(T)-D{i})/norm(D{i});
%               fprintf('tensor #%d has error %f\n',i,err);
%         >>end
%         tensor #2 has error 0.000000
%         tensor #3 has error 0.000000
%         tensor #4 has error 0.000000
%         tensor #5 has error 0.000000
%         tensor #6 has error 0.052674
%         tensor #7 has error 0.052758
%         tensor #8 has error 0.052781
%         tensor #9 has error 0.052790
%         tensor #10 has error 0.052791
%
%         % compute STA with forgetting factor alpha = 0.1
%         % very roughly speaking, 1/alpha is the equivalent window size,
%         % in this case, it is 10
%         >>[T,S] = STA(D{1},[2,3,4]);
%         >>alpha = 0.9;
%         >>for i = 2:10
%               [T,S] = STA(D{i},[2,3,4],T,S, alpha);
%               err = norm(full(T)-D{i})/norm(D{i});
%               fprintf('tensor #%d has error %f\n',i,err);
%         >>end
% 
%         tensor #2 has error 0.000000
%         tensor #3 has error 0.000000
%         tensor #4 has error 0.000000
%         tensor #5 has error 0.000000
%         tensor #6 has error 0.024951 % error increases when B starts to
%         tensor #7 has error 0.022921 % appear, but drops quickly because
%         tensor #8 has error 0.005310 % of the forgetting factorx
%         tensor #9 has error 0.005490
%         tensor #10 has error 0.001423
% 
%----------------------------------------------------------------
% Copyright: 2006
%         Jimeng Sun, Christos Faloutsos
%         All rights reserved.
% Please address questions or comments to: jimeng@cs.cmu.edu
%----------------------------------------------------------------
%% Input processing
if nargin<2
    error('not enough inputs');
end
N = ndims(Xnew);
if nargin<3 
    % no pior co-variance matrices, initialize to be truncated identity
    dv = size(Xnew);
    for i=1:N
        U{i} = eye(dv(i),R(i));
    end
else
    %check if the input R match the size of T.U
    for i=1:N
       if size(T.U{i},2)~=R(i)
           error('Rank mismatch on %d mode\n',i);
       end
    end
    U = T.U;
end
if nargin<4
    S = {};
    for i = 1:N
        S{i} = 1e-10*ones(R(i)); %no prior on energy; 
    end
end
if nargin<5
    alpha = 1; %no forgetting factor
end
if nargin<6
    samplingPercent = 0.1; %default sampling percentage is 10%
end
%% Doing STA for each mode
for n = 1:N
% matricize Xnew
    switch class(Xnew)
        case {'tensor'} %dense
            XM = double(tenmat(Xnew,n));
        case {'sptensor'} %sparse
            XM = double(sptenmat(Xnew,n));
    end
% Sampling the vectors
    ids = sampleCol(XM, samplingPercent);
% Update the projection matrices in the old tucker tensor
    Schange = zeros(size(S{n}));
    for i = ids
        x = XM(:,i);
        [U{n},S{n}] = STA_vec(x, S{n}, U{n}, alpha);        
    end
    for i = 1:size(U{n},2) %renormalization (precaution, maybe not necessary)
        U{n}(:,i) = U{n}(:,i)/norm(U{n}(:,i));
    end
end
%% compute core
core = ttm(Xnew, U, 't');
Tnew = ttensor(core, U);
return;
%% sampling function
function ids = sampleCol(X, p)
if nargin<2
    p = 1;
end
Xs = full(sum(abs(X),1));
[v, ids] = sort(Xs,'descend');
id2 = find(v>0);
id2 = id2(end);
ids = ids(1:max(floor(p*id2),1));
return;
%% STA update on a vector
function [U, dnew] = STA_vec(x, d, U, alpha)
%streamingTPCA: update U matrix based on vector x
    dnew = d;
    for i = 1:size(U,2) %update one column of U at a time
        [U(:,i),dnew(i), x] = updateW(x, U(:,i), d(i), alpha);
    end
    U = grams(U); % re-orthogonalization
return;
