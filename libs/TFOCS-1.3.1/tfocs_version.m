function tfocs_ver2 = tfocs_version

%TFOCS_VERSION   Version information.
%    Prints information about the version of TFOCS being used, and the
%    system on which it is running. When submitting a bug report, please
%    run this function and include its output in your report.

tfocs_ver = '1.3';
if nargout == 0,
    fprintf( 'TFOCS v%s\n', tfocs_ver );
    verd = ver('MATLAB');
    fprintf( 'MATLAB version %s %s on %s\n', verd.Version, verd.Release, computer );
else
    tfocs_ver2 = tfocs_ver;
end

% TFOCS v1.3 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2013 California Institute of Technology and CVX Research.
% See the file LICENSE for full license information.
