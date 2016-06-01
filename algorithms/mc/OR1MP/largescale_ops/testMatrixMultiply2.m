% test the sparse updates of augmented Lagrange.
function [] = testMatrixMultiply2()
clc;

data  = [2 3  4  5  6 ];
Known = [1 11 13 20 44];

m = 5; n = 10; r = 3;

[Ik, Jk] = ind2sub([m n],Known);
S_data   = sparse(Ik, Jk, data, m, n);
rho = 0.5;



% initially S_comp = empty. X0 = U0*V0 + S_comp
S_comp =  sparse(Ik, Jk, eps * ones(size(data)), m, n);

rng(100);
U0 = rand(m, r); U0s = U0;
V0 = rand(r, n); V0s = V0;

check_func = @(NameStr, X, Xs) fprintf('[%s] Related Error %.4g\n', NameStr, norm(X-Xs, 'fro')/norm(X, 'fro'));
%check_func = @(NameStr, X, Xs) fprintf('[%s] Related Error %.4g\n', NameStr, norm(X-Xs, 'fro'));

X0 = (U0 * V0 + S_comp);

disp('init lambda');
[X1  U1  V1]         = update_initLambda (X0, U0, V0, Known, data, rho);
[X1s U1s V1s S_comp] = update_initLambdaSparse (S_comp, U0s, V0s, Ik, Jk, data, rho, true);

check_func('U+', U1, U1s);
check_func('V+', V1, V1s);
check_func('X+', X1, X1s);
for i = 1: 2   
    [X1  U1  V1]         = update_initLambda (X1, U1, V1, Known, data, rho);
    [X1s U1s V1s S_comp] = update_initLambdaSparse (S_comp, U1s, V1s, Ik, Jk, data, rho, false);
    check_func('U+', U1, U1s);
    check_func('V+', V1, V1s);
    check_func('X+', X1, X1s);
end

disp('given sparse lambda');
U0 = rand(m, r); U0s = U0;
V0 = rand(r, n); V0s = V0;
X0 = (U0 * V0 + S_comp);

Lambda = sparse(Ik, Jk, ones(size(data)), m, n);
[X1 U1 V1] = update_givenLambda (X0, U0, V0, Known, data, rho, Lambda);
[X1s U1s V1s S_comp] = update_givenLambdaSparse (S_comp, U0s, V0s, Ik, Jk, data, rho, Lambda, true);
check_func('U+', U1, U1s);
check_func('V+', V1, V1s);
check_func('X+', X1, X1s);
for i = 1: 2   
    [X1  U1  V1]         = update_givenLambda (X1, U1, V1, Known, data, rho, Lambda);
    [X1s U1s V1s S_comp] = update_givenLambdaSparse (S_comp, U1s, V1s, Ik, Jk, data, rho, Lambda, false);
    check_func('U+', U1, U1s);
    check_func('V+', V1, V1s);
    check_func('X+', X1, X1s);
end



    function [X1 U1 V1] = update_initLambda (X0, U0, V0, Known, data, rho)
        Lambda0 = ones(m, n);
        % X, V => U
        T = X0 - Lambda0/rho;
        U1  = T * V0';
        % X, U => V
        V1  = U1' * T;        
        % U, V => X
        T2 = Lambda0/rho + U1 * V1;
        T3 = (2 * data + rho * T2(Known))/(2+ rho);
        X1 = T2; X1(Known) = T3;
    end

    function [X1s U1s V1s S_comp] = update_initLambdaSparse (S_comp, U0, V0, Ik, Jk, data, rho, first_iter)
        % Lambda is given by Lambda = e * et
        e  = ones(m, 1);
        et = ones(1, n);
        if first_iter
            % X, V => U
            U1s =  U0 * (V0 * V0') + S_comp  * V0' -  e * (et * V0') /rho;
            % X, U => V
            V1s = (U1s' * U0) * V0 + U1s' * S_comp  - (U1s' * e) * et / rho;
            % U, V => X
            s_comp_val = 2/(2+rho) * ( data - ones(size(data))/ rho - sparse_inp(U1s', V1s, Ik, Jk));
        else
            % X, V => U
            U1s =  U0 * (V0 * V0') + S_comp  * V0';
            % X, U => V
            V1s = (U1s' * U0) * V0 + U1s' * S_comp;
            % U, V => X
            s_comp_val = 2/(2+rho) * ( data  - sparse_inp(U1s', V1s, Ik, Jk));
        end
        sparse_update(S_comp, s_comp_val);
        X1s = e * et/rho + U1s * V1s + S_comp;
    end


    function [X1 U1 V1] = update_givenLambda (X0, U0, V0, Known, data, rho, Lambda)
        % X, V => U
        T = X0 - Lambda/rho;
        U1  = T * V0';
        % X, U => V
        V1  = U1' * T;        
        % U, V => X
        T2 = Lambda/rho + U1 * V1;
        T3 = (2 * data + rho * T2(Known))/(2+ rho);
        X1 = T2; X1(Known) = T3;
    end

    function [X1s U1s V1s S_comp] = update_givenLambdaSparse (S_comp, U0, V0, Ik, Jk, data, rho, LamSMat, first_iter)
        % Lambda is a sparse matrix here. 
        if first_iter
            % X, V => U
            U1s =  U0 * (V0 * V0') + S_comp  * V0' -  (LamSMat * V0') /rho;
            % X, U => V
            V1s = (U1s' * U0) * V0 + U1s' * S_comp  - (U1s' * LamSMat) / rho;
            % U, V => X
            [~, ~, ll] = find(LamSMat);
            s_comp_val = 2/(2+rho) * ( data - ll'/ rho - sparse_inp(U1s', V1s, Ik, Jk));
        else
            % X, V => U
            U1s =  U0 * (V0 * V0') + S_comp  * V0';
            % X, U => V
            V1s = (U1s' * U0) * V0 + U1s' * S_comp;
            % U, V => X
            s_comp_val = 2/(2+rho) * ( data  - sparse_inp(U1s', V1s, Ik, Jk));
        end
        sparse_update(S_comp, s_comp_val);
        X1s = LamSMat/rho + U1s * V1s + S_comp;
    end
end