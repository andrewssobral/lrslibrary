x = [9 .117512014694031425134547; ...
			2.5 .4903577561002348649728011; ...
			0.1 101.4332991507927588172155; ...
			7e-4 2040817.969783389022409943; ...
			7e-5 204081634.2978270192803090; ...
			7e-6 20408163266.95103968719027; ...
			7e-7 2040816326532.257177281929; ...
			-0.5 8.934802200544679309417246; ...
			-1.1 102.7489862404689390536693 ...
			];
for i = 1:rows(x)
	actual = trigamma(x(i,1));
	expected = x(i,2);
	e = abs(actual - expected)/expected;
	if e > 1e-12
		error(['trigamma(' x(i,1) ') = ' actual ' should be ' expected])
	end
end
if trigamma(-1) ~= -Inf
	error('trigamma(-1) should be -Inf');
end
if trigamma(0) ~= -Inf
  error('trigamma(0) should be -Inf');
end
if ~isnan(trigamma(-Inf))
  error('trigamma(-Inf) should be NaN');
end
