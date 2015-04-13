if 0
	s = sparse([1 2],[1 2],[2 3]);
	%awf_sparse(int32([1 2]),int32([1 2]),[2 3])
	setnonzeros(s,[4 5])
	return
end

n = 1000;
s = rand(n,n);
[i,j,v] = find(s);
[m,n] = size(s);
tic, for iter = 1:10, s = sparse(i,j,v,m,n); end; t1=toc;
fprintf('time for sparse = %g\n', t1);
%tic, for iter = 1:10, s = awf_sparse(int32(i),int32(j),v); end; t1=toc;
%fprintf('time for awf_sparse = %g\n', t1);
tic, for iter = 1:10, s = setnonzeros(s,v); end; t2=toc;
fprintf('time for setnonzeros = %g (%g times faster)\n', t2, t1/t2);
