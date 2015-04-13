function [compiler,options] = mexcompiler
% mexcompiler returns the name of the compiler used by mex.

% Written by Tom Minka

compiler = '';
options = struct;
options.vcvarsopts = '';
if ispc
	mexopts = fullfile(prefdir,'mexopts.bat');
	if ~exist(mexopts,'file')
		return
	end
	fid = fopen(mexopts);
	while 1
		txt = fgetl(fid);
		if ~ischar(txt), break, end
		if isempty(compiler) && strncmp('rem ',txt,4)
			compiler = txt(5:end-8);
		elseif ~isempty(strmatch('set ',txt))
			txt = txt(5:end);
			pos = strfind(txt,'=');
			if ~isempty(pos)
				pos = pos(1);
				field = txt(1:(pos-1));
				txt = txt((pos+1):end);
				txt = strrep(txt,'%MATLAB%',matlabroot);			
				options.(field) = txt;
				if strcmp(field,'VSINSTALLDIR')
					vsinstalldir = txt;
					if ~isempty(strmatch(':\Program Files (x86)\',vsinstalldir(2:end)))
						% http://msdn.microsoft.com/en-us/library/x4d2c09s(VS.80).aspx
						if exist(fullfile(vsinstalldir,'VC','bin','x86_amd64','cl.exe'))
							options.vcvarsopts = 'x86_amd64';
						else
							options.vcvarsopts = 'x86';
						end
					end
				end
			end
		end
	end
	fclose(fid);
else
	mexopts = fullfile(prefdir,'mexopts.sh');
	if ~exist(mexopts,'file')
		return
	end
	fid = fopen(mexopts);
	while 1
		txt = fgetl(fid);
		if ~ischar(txt), break, end
		if isempty(compiler) && strncmp('rem ',txt,4)
			compiler = txt(5:end-8);
		else
			pos = strfind(txt,'=');
			if ~isempty(pos)
				pos = pos(1);
				field = txt(1:(pos-1));
				txt = txt((pos+1):end);
				txt = strrep(txt,'%MATLAB%',matlabroot);			
				options.(field) = txt;
			end
		end
	end
	fclose(fid);
end
