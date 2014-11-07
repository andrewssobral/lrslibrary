%% projected dogleg method

g=proj_fr(vm.g(A,B),U,V);
g_norm=norm(g(:));

switch riop.sig
    case 1          % Riemannian gradient + Armijo linear search
        d=-g;
        [S1,U1,V1]=rtr_fr(d,S,U,V);
        A1=U1*S1*V1';
        
        v1=A_old(:)-A1(:);
        stat.ang1(end+1)=(vm.f(A_old,B)-vm.f(A1,B))/norm(v1(:))^2;
        
        if stat.ang1(end)<sg.thr
            for j=1:ls.imax
                d=d*.5;
                [S1,U1,V1]=rtr_fr(d,S,U,V);
                A1=U1*S1*V1';
                v1=A_old(:)-A1(:);
                stat.ang1(end)=(vm.f(A_old,B)-vm.f(A1,B))/norm(v1(:))^2;

                if stat.ang1(end)>=sg.thr
%                     disp([k j]);
                    break;
                end
            end
%             if j==ls.imax, disp(['linear search failed at iter ',num2str(k)]); end
        end
        
        S=S1; U=U1; V=V1; A=A1;
        
        
    case 2          % projected dogleg
        % Riemannian Hessian
%         H=@(d)reshape(riem_hess(reshape(d,n1,n2),vm.g(A,B),U,V,S),n1*n2,1);
        H=@(d)reshape(riem_hess(reshape(d,n1,n2),B-Z,U,V,S),n1*n2,1);
        
        % prepare dogleg
        gHg=g(:)'*H(g(:))/g_norm^2;
        if gHg<=riop.eps_c              % Cauchy point failed
            dogleg_sig=3;
            d=-g;
%             disp(['Cauchy point failed at iter ',num2str(k)])
        else
            d_c=-g/gHg;                 % Cauchy point
            [d_n,flag,relres,iter]=pcg(H,-g(:),cg.tol,cg.imax);
            stat.cg_count(end+1)=iter;
%             disp([k flag relres iter])
            d_n=reshape(d_n,n1,n2);     % Newton point
            if flag==4 || d_c(:)'*(d_n(:)-d_c(:))<0 %|| d_n(:)'*g(:)>=0  % negative curvature
                dogleg_sig=2;
                d=d_c;
%                 disp(['Negative curvature occurs at iter ',num2str(k)])
            else                % utilize dogleg method
                d=d_n;
                dogleg_sig=1;
            end
        end
%         disp([k (gHg-g_norm^2)/gHg])
        
        % dogleg search
        [S1,U1,V1]=rtr_fr(d,S,U,V);
        A1=U1*S1*V1';
        v1=A_old(:)-A1(:);
        stat.ang1(end+1)=(vm.f(A_old,B)-vm.f(A1,B))/norm(v1(:))^2;
        
        if stat.ang1(end)<sg.thr
            for j=1:ls.imax
                switch dogleg_sig
                    case 1
                        if j<=2
                            d=.5*j*d_c+(1-.5*j)*d_n;
                        else
                            d=d*.5;
                        end
                    case {2,3}
                        d=d*.5;
                end

                [S1,U1,V1]=rtr_fr(d,S,U,V);
                A1=U1*S1*V1';
                v1=A_old(:)-A1(:);
                stat.ang1(end)=(vm.f(A_old,B)-vm.f(A1,B))/norm(v1(:))^2;

                if stat.ang1(end)>=sg.thr
%                     disp([k j]);
                    break;
                end
            end
            stat.ls_count(end+1)=1+j;
%             if j==ls.imax, disp(['linear search failed at iter ',num2str(k)]); end
        else
            stat.ls_count(end+1)=1;
        end
        
        S=S1; U=U1; V=V1; A=A1;
end
