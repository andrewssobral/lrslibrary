function args = get_varargin(vars)
args.foo = 0;
for i=1:2:size(vars,2)
  args = setfield(args,vars{i},vars{i+1});
end

