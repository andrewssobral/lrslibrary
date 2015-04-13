% Add Manopt to the path and make all manopt components available.

% This file is part of Manopt: www.manopt.org.
% Original author: Nicolas Boumal, Jan. 3, 2013.
% Contributors: 
% Change log: 
%   Aug.  7, 2013 (NB): Changed to work without the import command
%                      (new structure of the toolbox).
%   Aug.  8, 2013 (NB): Changed to use addpath_recursive, home brewed.
%   Aug. 22, 2013 (NB): Using genpath instead of homecooked
%                       addpath_recursive.

addpath(pwd);

% Recursively add Manopt directories to the Matlab path.
cd manopt;
addpath(genpath(pwd));
cd ..;
