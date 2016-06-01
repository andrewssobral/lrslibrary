function GCO_BuildLib(Options)
% GCO_BuildLib    Attempt to compile and link the GCO_MATLAB library.
%    GCO_BuildLib is used internally by all other GCO_MATLAB commands 
%    to recompile the wrapper library if it is not yet built. 
%
%    YOU DO NOT NEED TO EXPLICITLY CALL THIS FUNCTION, unless you want to
%    customise the build settings via GCO_BuildLib(Options).
%    Default options:
%      Options.Debug=0                % optimised, detailed checking disabled
%      Options.EnergyType='long long' % int64 energy accumulator C type
%      Options.EnergyTermType='int'   % int32 energy term C type; if not specified, same as EnergyType
%
%    Example:
%      % Enable detailed assertions (e.g. than energy does not go up
%      % during expansion) and use double-precision float energy terms.
%      GCO_BuildLib(struct('Debug',1,'EnergyType','double'));
%

if (nargin < 1)
    Options = struct();
end
if (~isfield(Options,'Debug')), Options.Debug = 0; end
if (~isfield(Options,'EnergyType')), Options.EnergyType = ''; end
if (~isfield(Options,'EnergyTermType')), Options.EnergyTermType = ''; end
if (~isfield(Options,'EnergyType32')), Options.EnergyType32 = 0; end
if (~isfield(Options,'Force')), Options.Force = 1; end

MEXFLAGS = '';
if (strcmp(computer(),'GLNXA64') || strcmp(computer(),'PCWIN64') || strcmp(computer(),'MACI64'))
    MEXFLAGS = [MEXFLAGS ' -largeArrayDims -DA64BITS'];
end
if (Options.Debug)
    MEXFLAGS = [MEXFLAGS ' -g'];
end
if (Options.EnergyType32)
    MEXFLAGS = [MEXFLAGS ' -DGCO_ENERGYTYPE32'];
end
if (Options.EnergyType)
    MEXFLAGS = [MEXFLAGS ' -DGCO_ENERGYTYPE=' Options.EnergyType];
end
if (Options.EnergyTermType)
    MEXFLAGS = [MEXFLAGS ' -DGCO_ENERGYTERMTYPE=' Options.EnergyTermType];
end
if (strcmp(computer(),'PCWIN')) % link with libut for user interruptibility
    MEXFLAGS = [MEXFLAGS ' -D_WIN32 "' matlabroot() '\extern\lib\win32\microsoft\libut.lib"' ];
elseif (strcmp(computer(),'PCWIN64'))
    MEXFLAGS = [MEXFLAGS ' -D_WIN64 "' matlabroot() '\extern\lib\win64\microsoft\libut.lib"' ];
else
    MEXFLAGS = [MEXFLAGS ' -lut' ];
end

LIB_NAME = 'gco_matlab';
GCOMATDIR = fileparts(mfilename('fullpath'));
GCODIR = fileparts(GCOMATDIR);
OUTDIR = [ GCOMATDIR filesep 'bin' ];
[status msg msgid] = mkdir(GCOMATDIR, 'bin'); % Create bin directory
addpath(OUTDIR);                              % and add it to search path
if (~Options.Force && exist('gco_matlab')==3)
    return;
end
clear gco_matlab;

mexcmd = ['mex ' MEXFLAGS ' -outdir ''' OUTDIR ''' -output ' LIB_NAME ' ' ];

% Append all source file names to the MEX command string
SRCCPP = { 
    [GCOMATDIR filesep 'gco_matlab.cpp'],
    [GCODIR filesep 'GCoptimization.cpp'],
    [GCODIR filesep 'graph.cpp'],
    [GCODIR filesep 'maxflow.cpp'],
    [GCODIR filesep 'LinkedBlockList.cpp']
    };
for f=1:length(SRCCPP)
    mexcmd = [mexcmd ' ''' SRCCPP{f} ''' '];
end

eval(mexcmd);  % compile and link in one step

end
