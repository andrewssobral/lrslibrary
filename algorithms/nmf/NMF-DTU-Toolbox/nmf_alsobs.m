function [W,H]=nmf_alsobs(X,K,maxiter,speak)
%
% NMF using alternating least squares with obtimal brain surgeon.
%
% INPUT:
% X (N,M) : N (dimensionallity) x M (samples) non negative input matrix
% K       : Number of components
% maxiter : Maximum number of iterations to run
% speak   : prints iteration count and changes in connectivity matrix 
%           elements unless speak is 0
%
% OUTPUT:
% W       : N x K matrix
% H       : K x M matrix
%
% Kasper Winther Joergensen
% Informatics and Mathematical Modelling
% Technical University of Denmark
% kwj@imm.dtu.dk
% 2006/11/16

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% User adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

print_iter = 50; % iterations between print on screen and convergence test
obs_steps  = 15; % number of OBS steps to run before truncation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test for negative values in X
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if min(min(X)) < 0
    error('Input matrix elements can not be negative');
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialize random W and H
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[D,N]=size(X);
W=rand(D,K);
H=rand(K,N);

% use W*H to test for convergence
Xr_old = W*H;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Alternating least squares with 
% optimal brain surgeon iterations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for n=1:maxiter,
    W = ((pinv(H*H')*H)*X')';
    %%% OSB %%%
    count = 1;
    invHesW = (H*H')^-1;
    while count < obs_steps
        if min(min(W))<-eps
            dw = zeros(size(W));
            for i=1:size(W,1)
                ei = (W(i,:)<0);
                if sum(ei)
                    w_neg = W(i,ei);
                    h = invHesW(ei,ei);
                    la = h^-1*w_neg';
                    dw(i,:) = -invHesW(:,find(ei))*la;
                end
            end
            W  = W + dw;
        end
        count = count+1;
    end
    %%% OSB END %%%
    W = (W>0).*W; % truncate negative elements
    W = W./repmat(sum(W),D,1); % normalize columns to unit length

    %%%%%%%%%%%%%%%% 
    % H update
    %%%%%%%%%%%%%%%%
    H=(W*pinv(W'*W))'*X;
    %%% OSB %%%
    count = 1;
    invHesH = (W'*W)^-1;
    while count < obs_steps
        if min(min(H))<-eps
            dh = zeros(size(H));
            for i=1:size(H,2)
                ei = (H(:,i)<0);
                if sum(ei)
                    h_neg = H(ei,i);
                    h  = invHesH(ei,ei);
                    la = h^-1*h_neg;
                    dh(:,i) = -invHesH(:,find(ei))*la;
                end
            end
            H  = H + dh;
        end
        count = count+1;
    end
    %%% END OBS %%%
    H=H.*(H>0); % truncate negative elements
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % print to screen
    %%%%%%%%%%%%%%%%%%%%%%%
    if (rem(n,print_iter)==0) & speak,
        Xr = W*H;
        diff = sum(sum(abs(Xr_old-Xr)));
        Xr_old = Xr;
        eucl_dist = nmf_euclidean_dist(X,W*H);
        errorx = mean(mean(abs(X-W*H)))/mean(mean(X));
        disp(['Iter = ',int2str(n),...
            ', relative error = ',num2str(errorx),...
            ', diff = ', num2str(diff),...
            ', eucl dist ' num2str(eucl_dist)])
        if errorx < 10^(-5), break, end
    end
end