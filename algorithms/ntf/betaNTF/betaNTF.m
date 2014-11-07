function [W,H,Q, Vhat] = betaNTF(V,K,varargin)
%------------------------------------------------------------------
% simple beta-NTF implementation
%
%     Decomposes a tensor V of dimension FxTxI into a NTF model :
%           V(f,t,i) = \sum_k W(f,k)H(t,k)Q(i,k)
%
%     by minimizing a beta-divergence as a cost-functions.
%     particular cases include : 
%         * beta = 2 : Euclidean cost function
%         * beta = 1 : Kullback-Leibler costfunction
%         * beta = 0 : Itakura-Saito  cost function
%
%     It is possible to weight the cost function by a weight factor P of
%     the same size as V.
%
%
% inputs : 
%   * V        : FxTxI  non-negative data to approximate
%   * K        : number of NTF components
%
%   * Optional arguments include : 
%       - W : FxK  spectral basis
%       - H : TxK  temporal activations
%       - Q : IxK  loading factors
%       - P : FxTxI weighting matrix, permits to weight the cost function
%             to optimize element-wise
%       - fixedW  : do not learn W, default=0
%       - fixedH  : do not learn H, default=0
%       - fixedQ  : do not learn Q, default=0
%       - nIter   : number of iterations, default=100
%       - beta    : beta divergence considered, default=0 (Itakura-Saito)
%       - display : display plots during optimization (default=0)
% 
%
% outputs
%   * W,H,Q    : NTF model parameters
%   * Vhat     : reconstruction
%
% Example : 
%   [W,H,Q, Vhat] = betaNTF(V,'W',myW,'fixedW',1)
%--------------------------------------------------------------------------
%                                 Antoine Liutkus, Telecom ParisTech, 2012

%always same random seed for study
randn('seed', 0);

%Size of signals and number of NTF components
[F,T,I] = size(V);


%Parsing input and initialization
p = inputParser;
p.addParamValue('W', [], @(x)isnumeric(x));
p.addParamValue('H', [], @(x)isnumeric(x));
p.addParamValue('Q', [], @(x)isnumeric(x));
p.addParamValue('P', [], @(x)isnumeric(x));
p.addParamValue('fixedW', 0, @(x)isnumeric(x)||islogical(x));
p.addParamValue('fixedH', 0, @(x)isnumeric(x)||islogical(x));
p.addParamValue('fixedQ', 0, @(x)isnumeric(x)||islogical(x));
p.addParamValue('nIter', 70, @(x)(isnumeric(x)&&(x==round(x))));
p.addParamValue('beta', 0, @(x)isnumeric(x));
p.addParamValue('display', 0, @(x)isnumeric(x));

p.KeepUnmatched = true;
p.parse(varargin{:})

W         = p.Results.W;
H         = p.Results.H;
Q         = p.Results.Q;
P         = p.Results.P;
fixedW    = p.Results.fixedW;
fixedH    = p.Results.fixedH;
fixedQ    = p.Results.fixedQ;
nIter     = p.Results.nIter;
beta      = p.Results.beta;
display   = p.Results.display;



% Avoid zero values
%-----------------------------------
V = max(V,eps);

% initializing empty parameters
if isempty(W), W=.5*(1.5*abs(randn(F,K))+.5); end;
if isempty(H), H=.5*(1.5*abs(randn(T,K))+.5)*mean(V(:)); end;
if isempty(Q), Q=.5*(1.5*abs(randn(I,K))+.5); end;

% sources spectrogram model
Vhat = zeros(size(V));
for j = 1:I,  Vhat(:,:,j) = W * diag(Q(j,:)) * H'; end



%Main loop
%-----------------------------------
str=[];
disp('NTF model estimation');
for it = 1:nIter
    %display
    fprintf(repmat('\b',1,length(str)));
    str=sprintf('       Iteration %2d of total %2d', it, nIter);
    fprintf('%s', str);
    
    %update of W
    if ~fixedW
        Wnum = zeros(size(W));
        Wdenum = zeros(size(W));
        for j = 1:I
            if isempty(P)
                Wnum = Wnum + (Vhat(:,:,j).^(beta-2).*V(:,:,j))*H*diag(Q(j,:));
                Wdenum = Wdenum + (Vhat(:,:,j).^(beta-1))*H*diag(Q(j,:));
            else
                Wnum = Wnum + (P(:,:,j).*(Vhat(:,:,j).^(beta-2)).*V(:,:,j))*H*diag(Q(j,:));
                Wdenum = Wdenum + (P(:,:,j).*(Vhat(:,:,j).^(beta-1)))*H*diag(Q(j,:));
            end
        end
        W = W.*Wnum./Wdenum;
        
        %Update of models
        for j = 1:I,  Vhat(:,:,j) = W * diag(Q(j,:)) * H'; end
    end
    
    
    %update of H
    if ~fixedH
        Hnum = zeros(size(H))';
        Hdenum = zeros(size(H))';
        for j = 1:I
            if isempty(P)
                Hnum = Hnum + diag(Q(j,:))*W'*(Vhat(:,:,j).^(beta-2).*V(:,:,j));
                Hdenum = Hdenum + diag(Q(j,:))*W'*(Vhat(:,:,j).^(beta-1));
            else
                Hnum = Hnum + diag(Q(j,:))*W'*(P(:,:,j).*(Vhat(:,:,j).^(beta-2)).*V(:,:,j));
                Hdenum = Hdenum + diag(Q(j,:))*W'*(P(:,:,j).*(Vhat(:,:,j).^(beta-1)));
            end
        end
        H = H.*(Hnum'./Hdenum');
        
        %Update of models
        for j = 1:I,  Vhat(:,:,j) = W * diag(Q(j,:)) * H'; end
    end
    
    if ~fixedQ
        for j=1:I, 
            if isempty(P)
                Qnum = W'*(Vhat(:,:,j).^(beta-2).*V(:,:,j))*H;
                Qdenum =(W'*(Vhat(:,:,j).^(beta-1))*H );
            else
                Qnum = W'*(P(:,:,j).*(Vhat(:,:,j).^(beta-2)).*V(:,:,j))*H;
                Qdenum =(W'*(P(:,:,j).*(Vhat(:,:,j).^(beta-1)))*H );
            end
            Q(j,:) = diag(diag(Q(j,:)) .* Qnum./Qdenum);
        end

        %Update of models
        for j = 1:I,  Vhat(:,:,j) = W * diag(Q(j,:)) * H'; end
    end
   
    if display
        figure(1)
        clf
        for j = 1:I
            subplot(I,2,2*j-1)
            imagesc(log(V(:,:,j)))
            title(sprintf('original data, channel %d',j));
            xlabel('t');
            ylabel('f');
            grid on
            axis xy
            
            subplot(I,2,2*j)
            imagesc(log(Vhat(:,:,j)));
            title(sprintf('model, channel %d, K=%d components',j,K));
            xlabel('t');
            ylabel('f');
            grid on            
            axis xy
        end
        
        drawnow
    end
end

disp('   done.')
