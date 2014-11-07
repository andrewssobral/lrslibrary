function datadisp(T, dimlabels, opts)
%DATADISP Special display of a ktensor.
%
%   DATADISP(T,LABELS) displays the largest positive entries of each rank-1
%   factor of T using the corresponding labels. LABELS is a cell array of
%   size ndims(T) such that LABELS{n} is a string cell array of length
%   size(T,n). 
%
%   DATADISP(T,LABELS,OPTS) specify options:
%   OPTS.dimorder: Order to display the dimensions of T {1:ndims(T)}
%   OPTS.maxentries: Number of entries to show for each factor {10}
%   OPTS.printneg: Boolean to print the most negative entries {false}
%   OPTS.threshold: Threshold of smallest magnitude score to show {1e-4}
%
%   See also KTENSOR.
%
%MATLAB Tensor Toolbox.
%Copyright 2012, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.
% http://www.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2012) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in the file LICENSE.txt


%% Fill in optional variable
if ~exist('opts','var')
    opts = struct;
end

%% Set options from input or use defaults
dimorder = setparam(opts,'dimorder',1:ndims(T));
maxentries = setparam(opts,'maxentries',10);
printneg = setparam(opts,'printneg',false);
threshold = setparam(opts,'threshold',1e-6);

%% Main loop
R = size(T.lambda,1); % Rank
r = 1;
while (r <= R)

    fprintf(1, '\n======== Group %d ========\n', r);
    fprintf('\nWeight = %f\n', T.lambda(r));

    for i = dimorder(1:end)

        print_sublist(T.u{i}(:,r), dimlabels{i}, 'positive', maxentries, threshold);

        if printneg
            print_sublist(T.u{i}(:,r), dimlabels{i}, 'negative', maxentries, threshold);
        end

    end

    if r == R, 
        break, 
    end;

    foo = input('\nReturn to continue, jump to rank, or ''0'' (zero) to quit: ');
    if foo == 0
        return;
    elseif isempty(foo)
        r = r+1;
    else
        r = foo;
    end
end

return;


%%
function print_sublist(score, labels, type, maxentries, threshold)

if isequal(type,'positive')
    [sortedScore, sortedIdx] = sort(score, 'descend');
elseif isequal(type, 'negative')
    [sortedScore, sortedIdx] = sort(score, 'ascend');
else
    error('Invalid type');
end

sortedRefs = labels(sortedIdx);
entries = min([maxentries, length(score)]);

if isequal(type,'positive')
    range = 1:entries;
else
    range = entries:-1:1;
end

fprintf('%-10s %-4s %s\n','Score','Id','Name');

for k = range
    if abs(sortedScore(k)) < threshold
        continue;
    end
    if isequal(type,'negative') && sortedScore(k) >= 0
        continue;
    end
    if (abs(sortedScore(k)) < 1e-4)
        fprintf(1, '%10.3e  %4d %s\n', sortedScore(k), sortedIdx(k), ...
            sortedRefs{k});
    else
        fprintf(1, '%10.7f  %4d %s\n', sortedScore(k), sortedIdx(k), ...
            sortedRefs{k});
    end
end

%%
function x = setparam(opts,name,default)
if isfield(opts,name);
    x = opts.(name);
else
    x = default;
end
