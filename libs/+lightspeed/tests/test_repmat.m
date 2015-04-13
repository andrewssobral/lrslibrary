if 0
  % functionality test
  a = repmat([1 2; 3 4],2,4,2,2);
  b = repmat([1 2; 3 4],[2 4 2 2]);
	d = (a == b);
	assert(all(d(:)))
  repmat([1 2; 3 4],2)
  repmat(sqrt(-1),3,2)
  repmat('hello',3,2)
  repmat({7},3,2)
  repmat(sparse(7),3,2)
  s = struct('m',3);
  repmat(s,3,2)
  %x = Table(1,Universe(2));
  %repmat(x,3,2)
end

% run it once to load the definition
repmat(1,1,1);
if 0
  %x = rand(300,1);
  x = rand(10,1);
  niter = floor(1000000/prod(size(x)));
  n = 100;
  fprintf('repmat(rand(%g,%g),1,%g), %g times\n',rows(x),cols(x),n,niter);
  tic; for i = 1:niter a0=xrepmat(x,1,n); end; t0=toc; 
  fprintf('old repmat: %g\n',t0); 
  tic; for i = 1:niter a=repmat(x,1,n); end; t=toc; 
  assert(all(all(a0 == a)));
  fprintf('new repmat: %g (%g times faster)\n',t,t0/t);
end

if 0
  fprintf('old repmat:');
  tic; for i = 1:niter xrepmat(x',n,1); end; toc
  fprintf('new repmat:');
  tic; for i = 1:niter repmat(x',n,1); end; toc
end

if 0
  fprintf('old repmat:');
  tic; for i = 1:niter xrepmat(x,n,1); end; toc
  fprintf('new repmat:');
  tic; for i = 1:niter repmat(x,n,1); end; toc
end

if 0
  % repmat is faster than ones (in matlab v6)
  niter = 1000;
  fprintf('ones(300,1000):');
  tic; for i = 1:niter a0=ones(300,1000); end; toc
  fprintf('xones(300,1000):');
  tic; for i = 1:niter a=xones(300,1000); end; toc
  assert(all(all(a0==a)));
  fprintf('old repmat:');
  tic; for i = 1:niter a0=xrepmat(ones(300,1),1,1000); end; toc
  fprintf('new repmat:');
  tic; for i = 1:niter a=repmat(1,300,1000); end; toc
  assert(all(all(a0==a)));
end

if 0
  % zeros is faster than repmat (as expected)
  niter = 1000;
  fprintf('zeros(300,1000):');
  tic; for i = 1:niter zeros(300,1000); end; toc
  fprintf('old repmat:');
  tic; for i = 1:niter xrepmat(0,300,1000); end; toc
  fprintf('new repmat:');
  tic; for i = 1:niter repmat(zeros(300,1),1,1000); end; toc
end

if 1
  % test on large matrix
  b = rand(10240,1);
  n = 1;
	%n = 1;
	niter = 1000;
	iters = 1:niter;
  fprintf('repmat(rand(%g,%g),1,%g), %g times:\n',rows(b),cols(b),n,niter);
	if(1 || ~exist('t0')) 
		tic; for i=iters, a0=xrepmat(b,1,n); end; t0=toc;
		fprintf('old repmat: %g\n',t0); 
	end
  tic; for i=iters, a=repmat(b,1,n); end; t=toc;
  %assert(all(all(a0 == a)));
  fprintf('new repmat: %g (%g times faster)\n',t,t0/t);
	if n == 1
		tic; for i=iters, a=b+0; end; t1=toc;
		fprintf('plus: %g (%g times faster)\n',t1,t0/t1);
		tic; for i=iters, a=ones(10240,1); end; t2=toc;
		fprintf('ones: %g (%g times faster)\n',t2,t0/t2);
	end
end
if 0
	b = rand(10240,1);
	niter = 10240;
  tic; for i=1:niter, b+0; end; toc
  tic; for i=1:niter, a=b+0; end; toc
  tic; for i=1:niter, b+0; end; toc
  tic; for i=1:niter, a=b+0; end; toc
end
