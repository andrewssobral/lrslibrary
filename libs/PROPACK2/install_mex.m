function install_mex(VERB)
% install_mex
%   installs mex files for the PROPACK subpackage
% install_mex( verbose )
%   will display more verbose diagnoistic messages
%   if the variable "verbose" is true.
%
% Stephen Becker, Jan 2010, SVT package. srbecker@caltech.edu
%
%If you are having troubles with the compilation,
% make sure to run "mex -setup" to tell Matlab where
% your C compiler (e.g. gcc) is located. You only need to do this once.
%Matlab usually comes with a C compiler, but we know of the
%following exceptions:
%
%64-bit windows, Matlab versions 2008b and newer:
%    You will need to download an extern C/C++ compiler.
%    Information on a free C/C++ compiler is available here:
%    http://www.mathworks.com/support/solutions/en/data/1-6IJJ3L/index.html?solution=1-6IJJ3L
%
%64-bit Mac, recent Matlab versions:
%    You may need a C/C++ compiler.  You can get the gcc compiler
%    for free as part of Apple's "XCode" package
%        (you first need to signup for a free developer account):
%    http://developer.apple.com/xcode/)


if nargin < 1 || ~VERB
    VERBOSE = '';
else
    VERBOSE = '-v';
end

X='';
if ispc
    cc = get_compiler_config();
    if strcmpi(cc,'microsoft'),         cc = 'microsoft';
    else,         cc = 'lcc';     % Matlab's bundled compiler
    end
    c = computer;
    if strfind(c,'64')
        libpath = fullfile(matlabroot,'extern','lib','win64',cc);
    else
        libpath = fullfile(matlabroot,'extern','lib','win32',cc);
    end
    LAPACK = fullfile(libpath,'libmwlapack.lib');
    BLAS = fullfile(libpath,'libmwblas.lib');
    WIN = '-DWINDOWS';
else
    WIN = '-UWINDOWS';
    LAPACK = '-lmwlapack';
    % on linux/unix, sometimes MATLAB won't install its mwblas library
    % if the system already has a blas library.  So check for this:
    % (note, not in the same location as on Windows)
    c = lower(computer);
    if ismac, suffix = '.dylib'; else suffix = '.so'; end
    blasFile = fullfile(matlabroot,'bin',c, ['libmwblas',suffix] );
    if exist(blasFile,'file')
        BLAS = '-lmwblas';
    else
        BLAS = '-lblas';
        X='-DNO_BLAS';  % tell it not to include blas.h
    end
end
    

EXT = [];  % for native.  Use this when compiling for your own computer
% 2006a is v 7.2
if verLessThan('matlab', '7.3')
    LARGEARRAYDIMS = [];
    Y='-DNO_MATRIX_H'; % don't include matrix.h 'cause it doesn't exist!
else
    LARGEARRAYDIMS = '-largeArrayDims';
    Y=[];
end

OPT = '-O';

% compile PROPACK
mexHelper(VERBOSE,WIN,OPT,LARGEARRAYDIMS,EXT,'bdsqr_mex.c','dbdqr.c','-output','bdsqr',LAPACK,BLAS,X,Y);
mexHelper(VERBOSE,WIN,OPT,LARGEARRAYDIMS,EXT,'reorth_mex.c','reorth.c','-output','reorth',LAPACK,BLAS,X,Y);


function cc = get_compiler_config()
    % tested on Windows w/ R2008 only
    % This has to be in a function, otherwise old versions of matlab
    % get confused because "mex" is used as a structure (well, a class)
    % AND as a function.
    try 
        % this requires both a new verson of matlab and
        % that a compiler has been selected
        cc = mex.getCompilerConfigurations('C');
        cc = cc.Manufacturer;
        % Watch out for this error later (in old versions of matlab)
%         MATLAB:mir_error_function_previously_indexed_by_dot
    catch
        cc = [];
        fprintf('You may want to run ''mex -setup'' to setup the mex compiler,\n if you''ve never used the mex compiler before\n');
        wbsite='http://www.mathworks.com/support/solutions/en/data/1-6IJJ3L/index.html?solution=1-6IJJ3L';
        fprintf('If you have version 2008b or newer, and 64-bit Windows,\n');
        fprintf('then MATLAB does not come with a builtin compiler\n');
        fprintf('If you need a free C/C++ compiler, please see this mathworks website:\n%s\n',wbsite);
    end
