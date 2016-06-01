function p = get_params(params,name,defval)

if (isfield(params,name))
	p = getfield(params,name);
else
	p = defval;
end

