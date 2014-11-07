function [W, H, out] = LNMF(X, d, options)

% Spatially Localized Non-negative Matrix Factorization (LNMF)

% Reference:
%   S. Z. Li, X. Hou, H. Zhang, and Q. Cheng, "Learning spatially
%   localized, parts-based representation," in IEEE Conference on Computer Vision and Pattern Recognition, 2001, pp. 207–212.

% Written by Naiyang Guan (ny.guan@gmail.com)

[p, m] = size(X);

% Input arguments check out
if ~exist('options', 'var')
    options = [];
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

% avoid input error
X(X==0) = 1e-10;
X = X/max(X(:));

% test for negative values in X
if min(min(X)) < 0
    error('matrix entries can not be negative');
end
if min(sum(X,2)) == 0
    error('not all entries in a row can be zero');
end

niter = 3000;     % maximum number of iterations (can be adjusted)
precision = 1e-4;
f1 = 0;
nFreq = 10;

W = rand(p, d);
H = rand(d, m);
elapse = cputime;

for i=1:niter
    % divergence-reducing LNMF iterations
    W = (W./(repmat(sum(H, 2)', p, 1))).*((X./(W*H))*H');
    
    % normalization, suppose sigma(Bil)=1
%     B = B./(repmat(sum(B), p, 1));
    H = (H.*(W'*(X./(W*H)))).^0.5;
    
    % test convergence every 10 iterations
    if mod(i, nFreq)==1   
        % adjust small values to avoid undeflow
        W = max(W,eps);
        H=max(H,eps);
        f1 = KLC(X,W*H);
        if verbose,
            fprintf('%d,\t%.3f\n', i,f1);
        end
    end
    if mod(i, nFreq)==2
        f2 = KLC(X,W*H);
        if abs(f1-f2)<=precision
            break
        end
    end
end

out.ela = cputime - elapse;
out.itr = i;
out.fnm = norm(X-W*H, 'fro');
out.kld = KLC(X, W*H);

return;