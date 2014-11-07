function makeMex()
% Script to install mex files

% Change directory, so all mex files are in the mexFiles/ subdirectory
here = pwd;
cd( fullfile(tfocs_where,'mexFiles') );
if exist('OCTAVE_VERSION','builtin')
    octave = true; 
    compileFunction = @compileForOctave;
else
    octave = false;
    compileFunction = @compileForMatlab;
end


% -- Compile all the mex files --

% For the ordered L1 norm, used in prox_OL1.m:
compileFunction('proxAdaptiveL1Mex.c');




% Change directory back to wherever we started from
cd(here);

end




% -- subroutines --
function compileForMatlab( inputFile )
mex(inputFile)
end

function compileForOctave( inputFile )
if ~isunix,
    disp('warning, this will probably fail if not on linux; please edit makeMex.m yourself');
end
% To compile oct or mex files with octave, you must have the octave function "mkoctfile"
% e.g. (for ubuntu 10.04, using octave 3.2):   sudo apt-get install octave3.2-headers
%   for other versions of ubuntu and/or octave, try packages like
%       octave-pkg-dev or liboctave-dev
system(sprintf('mkoctfile --mex %s', inputFile ) );

% rm the .c ending
if length(inputFile) > 2 && inputFile(end-1:end)=='.c'
    inputFile = inputFile(1:end-2);
end
system(sprintf('rm %s.o', inputFile ) );
end


% TFOCS v1.3 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2013 California Institute of Technology and CVX Research.
% See the file LICENSE for full license information.
