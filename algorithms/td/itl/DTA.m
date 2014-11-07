function [T, Cnew] = DTA(Xnew, R,  C, alpha)
%Dynamic tensor decomposition DTA 
%   
%   [T, C] = DTA((Xnew, R, C, alpha) approximately updates tensor PCA,
%   according to new tensor Xnew, old variance matrices in C and he
%   specified dimensions in vector R. The input Xnew is a tensor,sptensor.
%   the result returned in T is a tucker tensor and C is the cell array of
%   new variance matrices.
%
%   Examples:
%         % construct two low rank tensor A and B
%         >>A = full(ttensor(tenrand([2,3,4]),{rand(10,2),rand(30,3), rand(40,4)})); 
%         >>B = full(ttensor(tenrand([2,3,4]),{rand(10,2),rand(30,3), rand(40,4)})); 
%         % create a tensor sequence/stream D with 5 copies of A and B
%         % respectively
%         >>D = {A,A,A,A,A,B,B,B,B,B};
%         % compute DTA with no forgetting factor
%         >>[T,C] = DTA(D{1},[2,3,4]);
%         >>for i = 2:10
%               [T,C] = DTA(D{i},[2,3,4],C);
%               err = norm(full(T)-D{i})/norm(D{i});
%               fprintf('tensor #%d has error %f\n',i,err);
%         >>end
%         tensor #2 has error 0.000000
%         tensor #3 has error 0.000000
%         tensor #4 has error 0.000000
%         tensor #5 has error 0.000000
%         tensor #6 has error 0.051471 % error increases when B 
%         tensor #7 has error 0.051435 % suddenly apprears instead of A
%         tensor #8 has error 0.051379
%         tensor #9 has error 0.051249
%         tensor #10 has error 0.050953
%
%         % compute DTA with forgetting factor alpha = 0.1
%         % very roughly speaking, 1/alpha is the equivalent window size,
%         % in this case, it is 10
%         >>[T,C] = DTA(D{1},[2,3,4]);
%         >>alpha = 0.1;
%         >>for i = 2:10           
%             [T,C] = DTA(D{i},[2,3,4],C, alpha);
%             err = norm(full(T)-D{i})/norm(D{i});
%             fprintf('tensor #%d has error %f\n',i,err);
%         >>end
%         tensor #2 has error 0.000000
%         tensor #3 has error 0.000000
%         tensor #4 has error 0.000000
%         tensor #5 has error 0.000000
%         tensor #6 has error 0.050496  % error increases when B starts to 
%         tensor #7 has error 0.049684  % appear, but quickly drops due to 
%         tensor #8 has error 0.008469  % the forgetting factor
%         tensor #9 has error 0.003901
%         tensor #10 has error 0.000078
%
%----------------------------------------------------------------
% Copyright: 2006
%         Jimeng Sun, Christos Faloutsos
%         All rights reserved.
% Please address questions or comments to: jimeng@cs.cmu.edu
%----------------------------------------------------------------
if nargin<2
    error('not enough inputs');
end
N = ndims(Xnew);
if nargin<3 % no pior co-variance matrices, initialize to be 0
    C = {};
    dv = size(Xnew);
    for n = 1:N
        C{end+1} = sparse(dv(n),dv(n));
    end
end
if nargin<4
    alpha = 1; %no forgetting factor
end
    
%% Update covariance matrices C and compute U
U = {};
for n = 1:N
    switch class(Xnew)
        case {'tensor'} %dense
            XM = double(tenmat(Xnew,n));
        case {'sptensor'} %sparse
            XM = double(sptenmat(Xnew,n));
    end
    Cnew{n} = alpha*C{n} + XM*XM';
    opts.disp = 0;
    opts.issym = 1;
    [U{n}, D] = eigs(Cnew{n},R(n),'LM',opts);    
end

%% compute core
core = ttm(Xnew, U, 't');
T = ttensor(core, U);
