function lambda = lambda_search(x, L, delta, xi)
% lambda = phi^{-1}(delta), 
% where phi(lambda)=|min{xi/lambda*E, (L/L+lambda)x}|_F and E=ones(m,n)

    treshold = 1e-12;
    norm_x = norm(x);
    if norm_x <= delta
        lambda = 0;
    else
        n = size(x,1);
        x = abs(x);
        x = sort(x, 'descend');
        
        if norm_x*(1-xi/(x(1)*L))<delta %This also includes the case: x(k)-xi/L<0
            lambda = L*(norm_x/delta - 1);
        else
            partial_sum = norm_x^2 - x(1)^2;
            flag = 0;
            for k=2:n
                xi_lambda_ratio = x(k)-xi/L;
                if xi_lambda_ratio > 0
                    tail_sum = (k-1)*(xi_lambda_ratio)^2;
                    test_sum = sqrt(tail_sum + (1-xi/(x(k)*L))^2*partial_sum);
                    if test_sum<delta
                        %lambda = L*(sqrt(partial_sum/(delta^2-tail_sum))-1);
                        coef = [delta^2, 2*L*delta^2, (delta^2-partial_sum)*L^2-(k-1)*xi^2, -2*L*(k-1)*xi^2, -(k-1)*xi^2*L^2];
                        lambda = roots(coef);
                        ind=(abs(imag(lambda))<treshold).*(lambda>0);
                        lambda = lambda'*ind;
                        flag = 1;
                        break
                    end
                    partial_sum = partial_sum - x(k)^2;
                else
                    coef = [delta^2, 2*L*delta^2, (delta^2-partial_sum)*L^2-(k-1)*xi^2, -2*L*(k-1)*xi^2, -(k-1)*xi^2*L^2];
                    lambda = roots(coef);
                    ind=(abs(imag(lambda))<treshold).*(lambda>0);
                    lambda = lambda'*ind;
                    flag = 1;
                    break
                end
            end
            if flag == 0
                lambda = sqrt(n)*xi/delta;
            end
        end
    end

    

    