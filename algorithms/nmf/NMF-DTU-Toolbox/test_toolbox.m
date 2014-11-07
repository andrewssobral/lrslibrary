clear

% create data
%X = rand(1000, 500);
load text_demo/exp_data_med;
X = sparse(X);

K = 5;
maxiter = 200;

% run nmf
alg = {'mm', 'cjlin', 'als', 'alsobs', 'prob'};
%alg = {'prob'}
for i=1:length(alg)
    [W{i} H{i}] = nmf(X, K, alg{i}, maxiter, 1);
end

% calc error
for i=1:length(alg)
   dist(i)=nmf_euclidean_dist(X,W{i}*H{i});
end
[y index] = sort(dist);

% display error
for i=index
   disp(['Euclidean distance ' num2str(dist(i),'%0.2f') ' for algorithm ', alg{i}]);
end
