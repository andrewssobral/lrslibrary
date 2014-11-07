%% alternating minimization scheme for low-rank+sparse decomposition

Z_norm=norm(Z(:));
stat.ener=[]; stat.ener(1)=vm.f(A,B);
if input_sig==1
A_norm=norm(A_0(:));
B_norm=norm(B_0(:));
stat.err_a=[]; stat.err_a(1)=norm(A(:)-A_0(:))/A_norm;
stat.err_b=[]; stat.err_b(1)=norm(B(:)-B_0(:))/B_norm;
end
stat.f_a=[]; stat.f_a(1)=norm(proj_fr(A+B-Z,U,V),'fro');
% stat.f_b=[]; %stat.f_b(1)=norm(proj_l0(A+B-Z,B),'fro');
% stat.ts1=[]; stat.ts2=[];
% stat.vi1=[]; stat.vi2=[];
stat.ang1=[]; stat.ang2=[];
stat.ls_count=[];
stat.cg_count=[];
stat.n3=[]; stat.n3(1)=n3;
stat.n4=[]; stat.n4(1)=n4;
stat.sfg_a=[]; stat.sfg_b=[];
stat.cpu=[]; %stat.cpu=time_main;
if ~isempty(A_ref), stat.diff_a=[]; stat.diff_a=norm(A(:)-A_ref(:))/A_norm; end
if ~isempty(B_ref), stat.diff_b=[]; stat.diff_b=norm(B(:)-B_ref(:))/B_norm; end
A_old=A; B_old=B;

% I1=eye(n1); I2=eye(n2);
waiting=waitbar(0,'iterating...');
%time_main=tic;
for k=1:slv.kmax
    waiting=waitbar(k/slv.kmax);
    
    % A-subproblem 
    if sg.sig~=1
    % first shot via partial SVD
    [U,S,V]=lansvd(Z-B,n3);
    A=U*S*V';
    
    % safeguard by projected dogleg method
    v1=A_old(:)-A(:);
    stat.ang1(end+1)=(vm.f(A_old,B)-vm.f(A,B))/norm(v1(:))^2;
    
    if sg.sig==2 && stat.ang1(end)<sg.thr
        A=A_old;
        stat.sfg_a(end+1)=1;
        disp(['A safeguarded by Riemannian optim at iter ',num2str(k)])
        riem_optim
    else
        stat.sfg_a(end+1)=0;
    end
    
    else    % directly via projected dogleg method
        riem_optim
    end
    
    % B-subproblem by sorting
    B=rtr_l0(Z-A,n4);
    
    % safeguard
    v1=B_old(:)-B(:);
    stat.ang2(end+1)=(vm.f(A,B_old)-vm.f(A,B))/norm(v1(:))^2;
    
    if sg.sig~=0 && stat.ang2(end)<sg.thr
        B=Z-A; B(B_old==0)=0;
        stat.ang2(end)=.5;
        stat.sfg_b(end+1)=1;
        disp(['B safeguarded at iter ',num2str(k)])
    else
        stat.sfg_b(end+1)=0;
    end
    
%     fprintf('|A^k-A^{k-1}|=%1.4e\n',norm(A(:)-A_old(:)))
    
    
    % trimming
    trimming
    
    
    % stats
    %stat.cpu(end+1)=stat.cpu(1)+toc(time_main);
    stat.ener(end+1)=vm.f(A,B);
    if input_sig==1
    stat.err_a(end+1)=norm(A(:)-A_0(:))/A_norm;
    stat.err_b(end+1)=norm(B(:)-B_0(:))/B_norm;
    end
    stat.f_a(end+1)=norm(proj_fr(vm.g(A,B_old),U,V),'fro');
%     stat.f_b(end+1)=norm(proj_l0(A+B-Z,B),'fro');
    if ~isempty(A_ref), stat.diff_a(end+1)=norm(A(:)-A_ref(:))/A_norm; end
    if ~isempty(B_ref), stat.diff_b(end+1)=norm(B(:)-B_ref(:))/B_norm; end
    stat.n3(end+1)=n3; stat.n4(end+1)=n4;
    
%     if sg.sig==1, cg.tol=1e-2*(stat.f_a(end)/stat.f_a(1))^.2; end
    
    % stopping rule
    if stat.f_a(end)<slv.ktol*stat.f_a(1)
%     if stat.err_a(end)<slv.ktol
% 	if stat.ener(end-1)-stat.ener(end)<slv.ktol*stat.ener(end-1)
        break
    end
    
    A_old=A;
    B_old=B;
end
%toc(time_main)
close(waiting)

% [sum(stat.ls_count),sum(stat.cg_count)]

