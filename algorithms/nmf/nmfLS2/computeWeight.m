function [weight] = computeWeight()
% compute distribution based weight for distributional 
% representation

% Load data
load 'sparseform.txt';
V = spconvert(sparseform);
method = 2;
fname = 'kld_weight.mat';

% Load files
load 'score.txt';

% Load sparse data
% score = train_score;
[R_idx, C_idx] = find(V);

% Data description
nSample = size(score, 1);
nFeature = size(V, 2);

% Initialize
% For positive class
% 1 - common 0's
% 2 - common 1's
% 3 - not common
% For negative class
% 4 - common 0's
% 5 - common 1's
% 6 - not common
featCount = zeros(6, nFeature);

% Counting numbers
for i = 1:nSample
    if mod(i, 100) == 0
        fprintf('.');
    end
    idx_1 = i*2-1; idx_2 = i*2;
    idx_feature_1 = C_idx(logical(R_idx==idx_1));
    idx_feature_2 = C_idx(logical(R_idx==idx_2));
    % common 0's, common 1's, not common
    [idx_c_0, idx_c_1, idx_u] = findIndex(idx_feature_1,...
        idx_feature_2, nFeature);
    if score(i) == 1
        featCount(1, idx_c_0) = featCount(1, idx_c_0) + 1;
        featCount(2, idx_c_1) = featCount(2, idx_c_1) + 1;
        featCount(3, idx_u) = featCount(3, idx_u) + 1;
    elseif score(i) == 0
        featCount(4, idx_c_0) = featCount(4, idx_c_0) + 1;
        featCount(5, idx_c_1) = featCount(5, idx_c_1) + 1;
        featCount(6, idx_u) = featCount(6, idx_u) + 1;
    end
end
fprintf('\n');
% Compute weights
featCount = featCount([2,3,5,6], : );
if method == 1
    weight = computeFeatureWeight(featCount);
elseif method == 2
    weight = computeFeatureWeight2(featCount);
elseif method == 3
    weight = computeFeatureWeight3(featCount);
end
save(fname, 'weight');
end

function [idx_0, idx_1, idx_u] = findIndex(idx_f1, idx_f2, nFeature)
tmp = zeros(nFeature, 2);
tmp(idx_f1', 1) = 1;
tmp(idx_f2', 2) = 1; 
idx_f1 = tmp(:,1);
idx_f2 = tmp(:,2);
% Common zeros
idx_0 = logical( (idx_f1==idx_f2) & (idx_f1==0) & (idx_f2==0) );
% Common ones
idx_1 = logical( (idx_f1==idx_f2) & (idx_f1==1) & (idx_f2==1) );
% Not common feature
idx_u = logical( idx_f1~= idx_f2 );
end

function [weight] = computeFeatureWeight(Freq)
% KL(P_S(X)||P_D(X))

Freq = Freq + 0.05;
% Pattern
pattern = [1,1,0,0; 1,1,0,0; 0,0,1,1; 0,0,1,1];
Freq = Freq./(pattern*Freq);
% 
ratio = log((Freq([1,2], :)./Freq([3,4], :)) + eps);
weight = ratio.*Freq([1,2], :);
weight = sum(weight);
end

function [weight] = computeFeatureWeight2(Freq)
% KL(P_D(X)||P_S(X))

Freq = Freq + 0.05;
% Pattern
pattern = [1,1,0,0; 1,1,0,0; 0,0,1,1; 0,0,1,1];
Freq = Freq./(pattern*Freq);
% 
ratio = log((Freq([3,4], :)./Freq([1,2], :)) + eps);
weight = ratio.*Freq([3,4], :);
weight = sum(weight);
end

function [weight] = computeFeatureWeight3(Freq)
% KL(P_S(X)||P_D(X))

Freq = Freq + 0.001;
% Pattern
pattern = [1,1,0,0; 1,1,0,0; 0,0,1,1; 0,0,1,1];
Freq = Freq./(pattern*Freq);
% Average
Freq_avg = (Freq([1,2],:) + Freq([3,4],:))/2;
% 
ratio = log((Freq([1,2], :)./Freq_avg)+ eps);
weight1 = ratio.*Freq([1,2], :);
weight1 = sum(weight1);
% 
ratio = log((Freq([3,4],:)./Freq_avg) + eps);
weight2 = ratio.*Freq([3,4], :);
weight2 = sum(weight2);
% 
weight = (weight1 + weight2)/2;
end
