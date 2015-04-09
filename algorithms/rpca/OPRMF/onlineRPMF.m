function [U, V, L] = onlineRPMF(X, rk, lambdaU, lambdaV, tol, mask)
	%% Initialize essencial variables
  %L = [];
	[m n] = size(X);
    maxIter = 40;
    % Use how many data to initialize
	%startIndex = 30;
  startIndex = 1;
	U = randn(m, rk);
	V = randn(rk, startIndex);
	lambda = 1;
	eps = 1e-3;
	IS = sparse(eye(rk));
    A = cell(n, 1);
    B = cell(n, 1);
    TA = cell(n, 1);
    TB = cell(n, 1);
    forgetFactor = 0.98; 
    % The first several estimation could be inaccuate, so we lower down
    % their confidence
    confidence = 1e-3;
	%% Main Process
	%tic
	for j = startIndex : n
		Y = X(:, 1 : j);
		if j ~= 1
            V(:, j) = V(:, j - 1);
		end
		r = abs(Y - U * V);
        % As more data come, the confidence should be bigger until 1.
        confidence = min(confidence * 1.5, 1);
		c = 0;
		disp(j);
		while true
			c = c + 1;
			oldR = r;
			%% Update V
			r = abs(Y - U * V);
			r = (r < eps) * eps + (r > eps) .* r;
			r = sqrt(lambda) ./ r; 
 
			if j == startIndex
				s = 1;
			else
				s = j;
			end
			parfor i = s : j
				T = U' * diag(sparse(r(:, i)) .* mask(:, i));
				V(:, i) = (T * U + lambdaV * IS) \ (T * Y(:, i));
            end
			%% Update U
			r = abs(Y - U * V);

			r = r';
			r = (r < eps) * eps + (r > eps) .* r;
			r = confidence * sqrt(lambda) ./ r;
            
			if j == startIndex
				parfor i = 1 : m
					T = V * diag(sparse(r(:, i)) .* mask(i, 1 : startIndex)');
					A{i} = inv(T * V' + lambdaU * IS);
					B{i} = T * Y(i, :)';
					U(i, :) = A{i} * B{i};
				end
            else
                v = V(:, j);
                TA = A;
                TB = B;
                parfor i = 1 : m
                    temp = A{i} * v / forgetFactor;
                    if mask(i, j) == 0
                        U(i, :) = TA{i} * TB{i};  
                        continue;
                    else
                        TA{i} = A{i} / forgetFactor - r(j, i) * (temp * temp') / (1 + r(j, i) * v' * temp);
                        TB{i} = B{i} * forgetFactor + r(j, i) * Y(i, j) * v;
                    end
                    U(i, :) = TA{i} * TB{i};  
                end
			end
 			r = abs(Y - U * V);
            if j == startIndex
               if sum(abs(r(:) - oldR(:))) / sum(oldR(:)) < tol && c ~= 1 || c > maxIter
                   L = U * V; 
                   break;
               end
            elseif sum(abs(r(:, j) - oldR(:, j))) / sum(oldR(:, j)) < tol || c > maxIter
                A = TA;
                B = TB;
                L(:, j) = U * V(:, j);
                break;
			end
        end
        disp(c);
    end
	%toc
end