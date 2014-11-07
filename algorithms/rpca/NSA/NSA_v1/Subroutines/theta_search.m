function theta = theta_search(x, rho, delta, xi)

% This function computes theta = phi^{-1}(delta), 
% where phi(theta)=|min{xi/theta*E, rho/(rho+theta)x}|_F and E=ones(m,n)

    treshold = 1e-12;
    norm_x = norm(x);
    if norm_x <= delta
        theta = 0;
    else
        n = size(x,1);
        x = abs(x);
        x = sort(x, 'descend');
        
        if norm_x*(1-xi/(x(1)*rho))<delta %This also includes the case: x(k)-xi/rho<0
            theta = rho*(norm_x/delta - 1);
        else
            partial_sum = norm_x^2 - x(1)^2;
            flag = 0;
            for k=2:n
                xi_theta_ratio = x(k)-xi/rho;
                if xi_theta_ratio > 0
                    tail_sum = (k-1)*(xi_theta_ratio)^2;
                    test_sum = sqrt(tail_sum + (1-xi/(x(k)*rho))^2*partial_sum);
                    if test_sum<delta
                        quartic_coef = [delta^2, 2*rho*delta^2, (delta^2-partial_sum)*rho^2-(k-1)*xi^2, -2*rho*(k-1)*xi^2, -(k-1)*xi^2*rho^2];
                        theta = roots(quartic_coef);
                        ind=(abs(imag(theta))<treshold).*(theta>0);
                        theta = theta'*ind;
                        flag = 1;
                        break
                    end
                    partial_sum = partial_sum - x(k)^2;
                else
                    quartic_coef = [delta^2, 2*rho*delta^2, (delta^2-partial_sum)*rho^2-(k-1)*xi^2, -2*rho*(k-1)*xi^2, -(k-1)*xi^2*rho^2];
                    theta = roots(quartic_coef);
                    ind=(abs(imag(theta))<treshold).*(theta>0);
                    theta = theta'*ind;
                    flag = 1;
                    break
                end
            end
            if flag == 0
                theta = sqrt(n)*xi/delta;
            end
        end
    end