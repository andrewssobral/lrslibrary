%% trimming

if trim.sig~=0 && mod(k,trim.mod)==0 
    % trim A via k-means
    v1=log10(diag(S));
    if norm(A(:)-A_old(:))/norm(A_old(:))<trim.tol && n3>1
        [i1,c1]=kmeans(v1,2);
        if max(c1)-min(c1)>trim.thr1
            i1=find(i1==i1(1));
            n3=i1(end);
            U=U(:,1:n3);
            V=V(:,1:n3);
            S=S(1:n3,1:n3);
            A=U*S*V';
            disp(['A trimmed at iter ',num2str(k)])
            if trim.disp
                figure(121); hold off;
                plot(v1,'r.','markersize',12); hold on;
                plot(v1(i1),'.','markersize',12); disp(c1)
                xlabel('iteration','fontsize',16)
                ylabel('singular value','fontsize',16)
                set(gca,'FontSize',16,'LineWidth',1)
            end
        end
    end

    % trim B
    v1=sort(abs(B(B~=0)),'descend');
    if norm(B(:)-B_old(:))/norm(B_old(:))<trim.tol && n4>1
        switch trim.sig
            case 1      % via k-means
                [i1,c1]=kmeans(v1,2);
%                 disp(max(c1)/min(c1))
                if max(c1)>min(c1)*trim.thr2
                    i1=find(i1==i1(1));
                    n4=i1(end);
                    B=rtr_l0(B,n4);
                    disp(['B trimmed at iter ',num2str(k)])
                    if trim.disp
                        figure(122); hold off;
                        plot(v1,'r.','markersize',12); hold on;
                        plot(v1(i1),'.','markersize',12); disp(c1)
                        xlabel('iteration','fontsize',16)
                        ylabel('absolute entry value','fontsize',16)
                        set(gca,'FontSize',16,'LineWidth',1)
                    end
                end
            
        case 2          % via hard thresholding
            if trim.disp
                figure(122); hold off;
                plot(v1,'r.','markersize',12); hold on;
                plot(v1(v1>=trim.thr2),'.','markersize',12);
                xlabel('iteration','fontsize',16)
                ylabel('absolute entry value','fontsize',16)
                set(gca,'FontSize',16,'LineWidth',1)
            end
            n4_old=n4;
            B(abs(B)<trim.thr2)=0;
            n4=sum(B(:)~=0);
            if n4<n4_old, disp(['B trimmed at iter ',num2str(k)]); end
        end
    end
end
