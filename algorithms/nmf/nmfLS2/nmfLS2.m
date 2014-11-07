function [W, H, err_ratio, err_v] = nmfLS2(V, K, max_iter,...
        use_nndsvd, tol, step)
if nargin < 3
    max_iter = 600;
end
if nargin < 4
    use_nndsvd = 1;
end
if nargin < 5
    tol = 1e-6;
end
if nargin < 6
    step=50;
end
% --------------------------------------------------
% Initialization
fprintf('Latent dimension is %d\n', K);
[nSample, nFeature] = size(V);
if use_nndsvd
    fprintf('Initialize W and H with NNDSVDa.\n')
    [W, H] = nndSVD(V, K);
    V_mean = mean(mean(V));
    W = W + V_mean;
    H = H + V_mean;
else
    fprintf('Initialize W and H randomly.\n')
    W = rand(nSample, K);
    H = rand(K, nFeature);
end
[W, H] = normalizeH(W, H, K);
[Row_idx, Col_idx, Val] = find(V);
Valprime = updateVprime(W, H, Row_idx, Col_idx);
% ---------------------------------------------------
% Updating
err_v = norm(Val-Valprime);
err_v_0 = err_v;
for iter = 1:max_iter
    if mod(iter, step) == 0
        fprintf('Iteration: %d\n', iter);
    end
    % --------------------------------------------
    % Update H 
    num = W'*V;
    det = (W'*W)*H;
    H = H.*(num./(det + eps));
    % Update W
    num = V*H';
    det = W*(H*H');
    W = W.*(num./(det + eps));
    % Normalize H
    [W, H] = normalizeH(W, H, K);
    % --------------------------------------------
    % Check the conditions
    if mod(iter, step) == 0
        err_v_OLD = err_v;
        Valprime = updateVprime(W, H, Row_idx, Col_idx);
        err_v = norm(Val-Valprime);
        % err_v_norm = norm(Val-Valprime);
        err_ratio = ((err_v_OLD-err_v)/(err_v_0-err_v));
        fprintf('The KL error is %f and the ratio is %f\n',...
            err_v, err_ratio);
        if err_ratio < tol
            fprintf('Program terminated by the given threshold (err=%f)\n', err_ratio);
            break;
        end
    end
    if iter >= max_iter
        fprintf('The number of iterations exceeds the max_iter\n');
        break;
    end
end
end

function [W, H] = normalizeH(W, H, K)
fprintf('Normalizing H ...\n');
weight = zeros(1,K);
for k = 1:K
    norm_h = norm(H(k,:), 2);
    H(k,:) = H(k,:)/norm_h;
    weight(k) = norm_h;
end
W = bsxfun(@times, W, weight);
end

function [w] = updateW(k, Val, W, H, Valprime, Row_idx, Col_idx)
Val = Val./(Valprime+eps);
h = H(k,:)'; % nFeature x 1
h_ext = h(Col_idx);
Val = h_ext.*Val;
V = sparse(Row_idx, Col_idx, Val);
Num = sum(V, 2); 
w = W(:, k).*(Num/sum(H(k,:)));
end

function [h] = updateH(k, Val, W, H, Valprime, Row_idx, Col_idx)
Val = Val./(Valprime+eps);
w = W(:, k);
w_ext = w(Row_idx);
Val = Val.*w_ext;
V = sparse(Row_idx, Col_idx, Val);
Num = sum(V, 1);
h = H(k, :).*(Num/sum(w));
end

function [Val] = updateVprime(W, H, Row_idx, Col_idx)
N = length(Row_idx);
Val = zeros(N,1);
for n = 1:N
    i = Row_idx(n);
    j = Col_idx(n);
    Val(n) = W(i,:)*H(:,j);
end
end
