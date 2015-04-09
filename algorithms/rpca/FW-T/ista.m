function [ output ] = ista( par )
% ista for stably principal component pursuit problem
%
%
%%  Cun Mu and John Wright, Mar '14

M = par.M; % data matrix
[m,n] = size(M); d = min(m,n);
Omega = ones(m,n); % if not specified, full observation
obj_target = -inf;

if isfield(par, 'Omega') Omega = par.Omega;  end
if isfield(par,'iter') iter_no = par.iter; end
if isfield(par,'obj') obj_target = par.obj; end

lambda_1 = par.lambda_1;
lambda_2 = par.lambda_2;
display = par.display; 
lips = 2;

L = zeros(m,n); S = zeros(m,n); 
sv = round(d/10); 
iter = 1; 
hist = 0; 

while true
    
    temp_L = L-Omega.*(L+S-M)/lips;
    temp_S = S-Omega.*(L+S-M)/lips;
     
    S = max(temp_S - lambda_2/lips, 0);
    S = S+min(temp_S + lambda_2/lips, 0);

    if choosvd(n, sv) == 1
        [U s V] = lansvd(temp_L, sv, 'L');
    else
        [U s V] = svd(temp_L, 'econ');
    end
    diagS = diag(s);
    svp = size(find(diagS > lambda_1/lips),1);
    if svp < sv
        flag = 1;
        sv = min(svp + 1, d);
    else
        flag = 0;
        sv = min(svp + round(0.05*d), d);
    end
    
    newSingularValues = diagS(1:svp) - lambda_1/lips;
    newNuclearNorm = sum(newSingularValues);
    
    L = U(:, 1:svp) * diag(diagS(1:svp) - lambda_1/lips) * V(:, 1:svp)';
    
    
    obj_value = 0.5*norm(Omega.*(L+S-M),'fro')^2 + lambda_1*newNuclearNorm ...
        +lambda_2*norm(vec(S),1);
    hist(iter) = obj_value;
      
    fprintf('this is the %d th iteration; sv = %d; obj. value = %d \n',...
        iter, sv, obj_value);
    
    if flag
        if obj_value < obj_target
            break;
        end
    end


    iter = iter +1;
    
end       

output.L = L;
output.S = S;
output.hist = hist;
end

