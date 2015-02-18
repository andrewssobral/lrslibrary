function [ output ] = FW_T( par)

% this function implements Frank-Wolfe-thresholding method
% 0.5*norm(Omega.*(L+S-M), 'fro')^2 + lambda_1 ||L||_* + lambda_2 ||S||_1
%
%
%%  Cun Mu and John Wright, Mar '14

M = par.M; % data matrix
[m,n] = size(M);

epsilon = 10^-3; compare = 0; % by default (if no specification)
if isfield(par, 'epsilon') epsilon = par.epsilon; end
lambda_1 = par.lambda_1;
lambda_2 = par.lambda_2;
iter = par.iter;
display = par.display;
method = par.method; % power or exact
Omega = par.Omega; 
rho = sum(sum(Omega))/m/n;
%% FW-T method

% initialization
L = zeros(m,n); S = zeros(m,n);
t_1 = 0; t_2 = 0;

U_1 = 0.5*norm(M,'fro')^2/lambda_1; % initial rough guess
U_2 = 0.5*norm(M,'fro')^2/lambda_2; % initial rough guess

history = 0; % to store the values of each iteration
tic
%% full observation scenario
%if ~isfield(par, 'Omega')
if rho == 1
    fprintf('full observation. \n');
    history(1) = norm(M,'fro')^2/2;
    temp = L + S - M;  % gradient
    count = 0;
    for k = 1: iter
        
        if k>=2
            if abs(history(k)-history(k-1))/history(k-1)<epsilon 
                count = count+1;
            else
                count=0;
            end
            if count==5
                break;
                fprintf('--------------------------------\n')
                fprintf('total # of iter. used is %d. \n', k);
            end
        end
        
        fprintf('the current function value is %d \n', history(k));
        
        if mod(k,1) == 0
            fprintf('--------------------------------\n')
            fprintf('this is the %d th iter. \n', k);
        end
        
        
        %------------------linearization subproblem-----------------------%
        
        if  strcmp(method,'power') % approximate
            fprintf('power method is leveraged \n')
            [U ev V] =  power_method(temp, 5);
        else % exact
            [U,ev,V] = lansvd(temp,1,'l'); % top eigen
        end
        
        D_L = -U*V';
        
        if lambda_1 >= ev
            V_L = 0; V_t_1 = 0;
        else
            V_L = U_1*D_L; V_t_1 = U_1;
        end
        
        [mag ind] = max(vec(abs(temp)));
        j = floor((ind-1)/m)+1; i = mod(ind-1,m)+1;
        sign_ = sign(temp(i,j));
        D_S = zeros(m,n);
        D_S(i,j) = -sign_;
        
        if lambda_2 >= mag
            V_S = 0; V_t_2 = 0;
        else
            V_S = U_2*D_S; V_t_2 = U_2;
        end
                
        H = zeros(2,2);
        temp_1 = V_L-L; temp_2 = V_S-S; % temp = L + S - M;
        H(1,1) = norm(temp_1, 'fro')^2;
        H(2,2) = norm(temp_2, 'fro')^2;
        H(1,2) = sum(sum(temp_1.*temp_2));
        H(2,1) = H(1,2);
        f = zeros(1,2);
        f(1) = sum(sum(temp_1.*temp));
        f(2) = sum(sum(temp_2.*temp));
        f = f + [lambda_1*(V_t_1-t_1), lambda_2*(V_t_2-t_2)];
        
        % using QP solvers
        lb = zeros(2,1);
        ub = [1;1];
        options.Display = 'off';
        options.TolFun = 10^-5;
        x = quadprog(H,f,[],[],[],[],lb,ub,[],options);
        x = quadprog(H,f,[],[],[],[],lb,ub,[],options);

        
        
        %----------------------- update L and S --------------------------%
        alpha = x(1);
        beta = x(2);
        L = (1-alpha)*L + alpha*V_L; t_1 = t_1 + alpha*(V_t_1-t_1);
        S = (1-beta)*S + beta*V_S; t_2 = t_2 + beta*(V_t_2-t_2);
       
        %------------------------ thresholding----------------------------%
        
        temp_3 = M-L;
        
        S = max(temp_3 - lambda_2, 0);
        S = S + min(temp_3 + lambda_2, 0);
        
        t_2 = sum(sum(abs(S)));
        
        %----------- update U_1 and U_2 to a better esitmate -------------%
        temp = L + S - M;
        history(k+1) = norm(temp,'fro')^2/2+lambda_1*t_1+lambda_2*t_2;
        U_1 = min(history(k+1)/lambda_1,U_1);
        U_2 = min(history(k+1)/lambda_2,U_2);
        
    end
    
end

%% partial observation scenario
if rho<1
    fprintf('partial observation. \n');
    temp = Omega.*(L + S - M);  % gradient
    history(1) = norm(temp,'fro')^2/2+lambda_1*t_1+lambda_2*t_2;
    count = 0;
    
    for k = 1: iter
        
        if k>=2
            if abs(history(k)-history(k-1))*2/history(1)<epsilon
                count = count+1;
            else
                count=0;
            end
            if count==5
                break;
                fprintf('--------------------------------\n')
                fprintf('total # of iter. used is %d. \n', k);
            end
        end
        
        % current function value
        fprintf('the current function value is %d \n', history(k));
        fprintf('--------------------------------\n')
        fprintf('this is the %d th iter. \n', k);
        
        
        
        %------------------linearization subproblem-----------------------%
        
        if  strcmp(method,'power') % approximate
            [U ev V] =  power_method(temp, 5);
        else % exact
            [U,ev,V] = lansvd(temp,1,'l'); % top eigen
        end
        
        D_L = -U*V';
        
        if lambda_1 >= ev
            V_L = 0; V_t_1 = 0;
        else
            V_L = U_1*D_L; V_t_1 = U_1;
        end
        
        [mag ind] = max(vec(abs(temp)));
        j = floor((ind-1)/m)+1; i = mod(ind-1,m)+1;
        sign_ = sign(temp(i,j));
        D_S = zeros(m,n);
        D_S(i,j) = -sign_;
        
        if lambda_2 >= mag
            V_S = 0; V_t_2 = 0;
        else
            V_S = U_2*D_S; V_t_2 = U_2;
        end
        
        
        %--------------------- use QP (exact search) ---------------------%
        H = zeros(2,2);
        temp_1 = Omega.*(V_L-L); temp_2 = Omega.*(V_S-S); % temp = L + S - M;
        H(1,1) = norm(temp_1, 'fro')^2;
        H(2,2) = norm(temp_2, 'fro')^2;
        H(1,2) = sum(sum(temp_1.*temp_2));
        H(2,1) = H(1,2);
        f = zeros(1,2);
        f(1) = sum(sum(temp_1.*temp));
        f(2) = sum(sum(temp_2.*temp));
        f = f + [lambda_1*(V_t_1-t_1), lambda_2*(V_t_2-t_2)];
        
        % using QP solvers
        lb = zeros(2,1);
        ub = [1;1];
        options.Display = 'off';
        options.TolFun = 10^-5;
        x = quadprog(H,f,[],[],[],[],lb,ub,[],options);
        x = quadprog(H,f,[],[],[],[],lb,ub,[],options);
        
        
        
        %----------------------- update L and S --------------------------%
        alpha = x(1);
        beta = x(2);
        L = (1-alpha)*L + alpha*V_L; t_1 = t_1 + alpha*(V_t_1-t_1);
        S = (1-beta)*S + beta*V_S; t_2 = t_2 + beta*(V_t_2-t_2);
        
        %------------------------ thresholding----------------------------%
        temp_3 = S-Omega.*(L+S-M);;
        S = max(temp_3 - lambda_2, 0);
        S = S + min(temp_3 + lambda_2, 0);
        t_2 = norm(vec(S),1);
        
        %----------- update U_1 and U_2 to a better esitmate -------------%
        temp = Omega.*(L + S - M);  % gradient
        history(k+1) = norm(temp,'fro')^2/2+lambda_1*t_1+lambda_2*t_2;
        U_1 = min(history(k)/lambda_1,U_1);
        U_2 = min(history(k)/lambda_2,U_2);
        
    end
    
end
toc
%% output
output.L = L;
output.S = S;
output.hist = history;
output.iter = k;
end

