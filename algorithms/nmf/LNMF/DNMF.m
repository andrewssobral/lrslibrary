function [B, H, out] = DNMF(X, L, d, options)

% Discriminative Non-negative Matrix Factorization (DNMF)
%
% Reference:
%   S. Zafeiriou, A. Tefas, I. Buciu, and I. Pitas, "Exploiting
%   Discriminant Information in Nonnegative Matrix Factorization with Application to Frontal Face Verification,"
%   IEEE transactions on neural networks, vol. 17, no. 3, pp. 683–95, May
%   2006.
%
% Arguments
% <Input>:
%   X: data matrix (nSmp x nFea);
%   L: label of each observation;
%   d: lower dimensionality;
%   options:
%       gamma: trade-off of regularizations;
%       delta: trade-off between within-class and between-class scatters;
%   verbose: details,
% <Output>:
%   B: basis matrix;
%   H: coefficient matrix.
%
% Written by Naiyang Guan (ny.guan@gmail.com)

% data rescaling
X = max(eps,X);
X = X/max(X(:));
[p, m] = size(X);

% Input arguments check out
if ~exist('options', 'var')
    options = [];
end
if ~isfield(options, 'gamma') || ~isfield(options, 'delta')
    gma = 1;    dlt = 1e-3;
else
    gma = options.gamma;    dlt = options.delta;
end
if ~isfield(options, 'verbose') || ~options.verbose
    verbose = 0;
else
    verbose = 1;
    if ~isfield(options, 'len') || ~isfield(options, 'wid')
        len = sqrt(p);        wid = len;
    else
        len = options.len;        wid = options.wid;
    end
end

% test for negative values in X
if min(min(X)) < 0
    error('matrix entries can not be negative');
end
if min(sum(X,2)) == 0
    error('not all entries in a row can be zero');
end

stopconv = 40;      % stopping criterion (can be adjusted)
niter = 5000;       % maximum number of iterations (can be adjusted)

cons = zeros(m,m);
consold = cons;
inc = 0;

B = rand(p, d);
B = B./repmat(sum(B), p, 1);
H = rand(d, m);
A = CalculateA(H, L, gma, dlt);

discold = 0;
objold = 0;

elapse = cputime;
for i=1:niter
    % divergence-reducing DNMF iterations
    [T] = CalculateT(H, L, gma, dlt);
    H = (T + sqrt(T.^2 + 4*A.*H.*(B'*(X./(B*H)))))./(2*A);
    
    B = (B./(repmat(sum(H, 2)', p, 1))).*((X./(B*H))*H');
    % normalization
    B = B./(repmat(sum(B), p, 1));
    
    % test convergence every 10 iterations
    if(mod(i,10)==0)  
    
        % adjust small values to avoid undeflow
        B = max(B,eps); 
        H=max(H,eps);
        
        % construct connectivity matrix    
        [y, index] = max(H, [], 1);
        mat1 = repmat(index, m, 1);
        mat2 = repmat(index', 1, m);
        cons = (mat1==mat2);
    
        if(sum(sum(cons~=consold))==0)
            inc = inc + 1;
        else
            inc = 0;
        end
            % Draw objective function and print them out
        if verbose
            [Sw, Sb] = CalculateD(L,H);
            discold = [discold Sb/Sw];
            objold = [objold KL(X,B*H)+gma*Sw-dlt*Sb];
            % Plot out
            figure(1);
            subplot(2,2,1);
            plot(discold(2:end));    ylabel('tr(Sb)/tr(Sw)');
            subplot(2,2,3);
            plot(objold(2:end));     ylabel('objective');
            subplot(2,2,[2,4]);
            ShowImage(B, len, wid, 5, 5, 1);
            drawnow;
            
            fprintf('\t%d\t%d\t%d\t%.2f\t%.2f\t%.2f\n', i, inc, sum(sum(cons~=consold)), KL(X, B*H), Sw, Sb);
        end
    
        if(inc > stopconv)
            break;
        end
    
        consold = cons;
    end
end

out.ela = cputime - elapse;
out.itr = i;
out.fnm = norm(X-B*H, 'fro');
out.kld = KLC(X, B*H);

return;

function A = CalculateA(H, L, gma, dlt)

% construct A of the qurdratic equation for DNMF
% low-dimensional coordinates, H
% label of the high-dimensional samples, L
% auother: Guan Naiyang

A = ones(size(H));
label = unique(L);

% update class by class
for r = 1:length(label)
    [index] = find(L == label(r));
    S = H(:,index);
    A(:,index) = ones(size(S))/length(index);
end

% gather together
A = 2*gma - 2*(gma+dlt)*A;

return;

function [Sw, Sb] = CalculateD(L, H)

% Construct discriminant for DNMF and CSDNMF

S = zeros(size(H));    % between class scatter matrix
label = unique(L);
M = mean(H, 2);        % global mean

for r = 1:length(label)
    [index] = find(L == label(r));
    C = H(:, index);                             % subclass
    m = repmat(mean(C, 2), 1, length(index));    % local mean
    H(:, index) = C - m;
    S(:, index) = m - repmat(M, 1, length(index));
end

Sw = trace(H'*H);
Sb = trace(S'*S);

return;

function [T] = CalculateT(H, L, gma, dlt)

% construct T1~T3 for DNMF from label L
% low-dimensional coordinates, H
% label of the high-dimensional samples, L
% auother: Guan Naiyang

label = unique(L);
M = mean(H, 2);

% update class by class
for r = 1:length(label)
    [index] = find(L == label(r));
    S = H(:,index);
    H(:,index) = repmat(mean(S,2), 1, length(index)) - S/length(index);
end

% gather together
T = 2*(gma+dlt)*H - 2*dlt*repmat(M, 1, length(L)) - 1;

return;