function s = SetDefaultParams(s, p)
% s = SetDefaultParams(s);
% Sets default parameters
% s: user-specified parameters that are used instead of defaults
%Written by: Prateek Jain (pjain@cs.utexas.edu) and Raghu Meka (raghu@cs.utexas.edu)
%Last updated: October, 2009

if (isfield(s, 'tol') == 0),
    s.tol = 1e-2;
end

if (isfield(s, 'vtol') == 0),
    s.vtol = 1e-3;
end

if (isfield(s,'maxtol')==0),
    s.maxtol=1e+3;
end

if (isfield(s, 'mxitr') == 0),
    s.mxitr = 100;
end

if (isfield(s, 'verbosity') == 0),
    s.verbosity = 1;
end

if (isfield(s, 'eta') == 0),
    delta=1/3;
    s.eta=1/(1+delta)/p;
end

