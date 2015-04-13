pairs = [0.1 2.2527126517342059598697; ...
			0.6 .39823385806923489961685; ...
			0.7 .26086724653166651438573; ...
			1.0 0; ...
			2.0 0; ...
			3.4 1.0923280598027415674947; ...
			4.0 1.791759469228055000812477; ...
			8.0 8.525161361065414300165531; ...
			64.0 201.00931639928152667928; ...
			256.0 1161.71210111840065079];
err = [];
for i = 1:rows(pairs)
	err(i) = abs(gammaln(pairs(i,1)) - pairs(i,2))/pairs(i,1);	
end
if max(err) > 1e-10
	err
	error('maximum err > 1e-10')
end
err = [];
for i = 1:rows(pairs)
	err(i) = abs(gammaln(pairs(i,1),1) - pairs(i,2))/pairs(i,1);	
end
if max(err) > 1e-10
	err
	error('maximum err > 1e-10')
end
err = abs(gammaln(1.1,2) - 0.920726359734123);
if err > 1e-10
	error('gammaln(1.1,2) != 0.920726359734123');
end
if gammaln(0) ~= Inf
  error('gammaln(0) should be Inf');
end
if ~isnan(gammaln(-1))
  %error('gammaln(-1) should be NaN');
end
if gammaln(Inf) ~= Inf
  error('gammaln(Inf) should be Inf');
end
% should be NaN?
%gammaln(-Inf)
if ~isnan(gammaln(NaN))
  error('gammaln(NaN) should be NaN');
end

