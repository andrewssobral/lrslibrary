function [L,S,error,time]=GreGoDec(D,rank,tau,tol,power,k)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                 Greedy Semi-Soft GoDec Algotithm (GreBsmo)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUTS:
%X: nxp data matrix with n samples and p features
%rank: rank(L)<=rank
%tau: soft thresholding
%power: >=0, power scheme modification, increasing it lead to better
%k: rank stepsize
%accuracy and more time cost
%OUTPUTS:
%L:Low-rank part
%S:Sparse part
%RMSE: error
%error: ||X-L-S||/||X||
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%REFERENCE:
% Tianyi Zhou and Dacheng Tao, "GoDec: Randomized Lo-rank & Sparse Matrix
% Decomposition in Noisy Case", ICML 2011
% Tianyi Zhou and Dacheng Tao, "Greedy Bilateral Sketch, Completion and 
% Smoothing", AISTATS 2013.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Tianyi Zhou, 2013, All rights reserved.

%matrix size
[m,n]=size(D);
if m<n D=D'; end
normD=norm(D(:));

%initialization of L and S
rankk=round(rank/k);
error=zeros(rankk*power,1);
[X,s,Y]=svds(D,k,'L');
%[X,s,Y]=svds(D,k);
Y=Y';
X=X*s;
L=X*Y;
S=wthresh(D-L,'s',tau);
T=D-L-S;
error(1)=norm(T(:))/normD;
iii=1;
stop=false;
alf=0;

tic;
for r=1:rankk
    % parameters for alf
    rrank=rank;
    est_rank = 1;
    rank_min =  1;
    rk_jump = 10;
    alf=0;increment=1;
    itr_rank = 0; minitr_reduce_rank = 5;    maxitr_reduce_rank = 50;
    if iii==power*(r-2)+1
        iii=iii+power;
    end
    for iter=1:power
        %Update of X
        X=L*Y';
        if est_rank==1
            [X,R,E]=qr(X,0);
        else
            [X,R]=qr(X,0);
        end
        
        % Update of Y
        Y=X'*L;
        L=X*Y;
        
        %Update of S
        T=D-L;
        S=wthresh(T,'s',tau);
        
        %Error, stopping criteria
        T=T-S;
        ii=iii+iter;
        error(ii)=norm(T(:))/normD;
        if error(ii)<tol
            stop=true;
            break;
        end
        
        % adjust est_rank
        if est_rank >= 1; rank_estimator_adaptive(); end
        if rrank ~= rank; rank = rrank; if est_rank ==0; alf = 0; continue; end; end
        % adjust alf
        ratio=error(ii)/error(ii-1);
        if ratio >= 1.1
            increment = max(0.1*alf, 0.1*increment);
            X = X1; Y = Y1; L = L1; S = S1; T = T1; error(ii) = error(ii-1);
            alf = 0;
        elseif ratio > 0.7;
            increment = max(increment, 0.25*alf);
            alf = alf + increment;
        end
        
        % Update of L
        X1=X;Y1=Y;L1=L;S1=S;T1=T;
        L=L+(1+alf).*T;
        
        % Add coreset
        if iter>8
            if mean(error(ii-7:ii))/error(ii-8)>0.92
                iii=ii;
                sf=size(X,2);
                if size(Y,1)-sf>=k
                    Y=Y(1:sf,:);
                end
                break;
            end
        end
    end
    % Stop
    if stop
        break;
    end
    % Coreset
    if r<rankk
        %[u,s,v]=lansvd(L,k,'L');
        %Y=[Y;v'];
        %[u,s,v]=svds(L,k);
        %Y=[Y;v'];
                 v=randn(k,m)*L;
                 Y=[Y;v];
    end
end
time=toc;

error(error==0)=[];
L=X*Y;
if m<n 
    L=L';
    S=S';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function rank_estimator_adaptive()
        if est_rank == 1
            dR = abs(diag(R));       drops = dR(1:end-1)./dR(2:end);
            [dmx,imx] = max(drops);  rel_drp = (rank-1)*dmx/(sum(drops)-dmx);
            if (rel_drp > rk_jump && itr_rank > minitr_reduce_rank) ...
                    || itr_rank > maxitr_reduce_rank; %bar(drops); pause;
                rrank = max([imx, floor(0.1*rank), rank_min]);
                error(ii) = norm(res)/normz;
                est_rank = 0; itr_rank = 0;
            end
        end
    end
end