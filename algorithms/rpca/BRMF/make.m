function make
% detect platform 

compstr = computer;
is64bit = strcmp(compstr(end-1:end),'64');


% compilation parameters

compile_params = cell(0);
if (is64bit)
  compile_params{1} = '-largeArrayDims';
end


% Compile files %

sources = {'mexutils.c'};

cd mex
disp('Compiling sampleTau...');
mex('sampleTau.cpp', sources{:},compile_params{:});
cd ..