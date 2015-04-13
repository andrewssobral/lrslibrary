x = mutable;
x(1) = 5;
y = x;
x(1) = 4;
if y(1) ~= 4
	error('y(1) ~= 4')
end

x = mutable(struct);
y = x;
x.a = 4;
if y.a ~= 4
	error('y.a ~= 4')
end
fprintf('Test passed.\n');
