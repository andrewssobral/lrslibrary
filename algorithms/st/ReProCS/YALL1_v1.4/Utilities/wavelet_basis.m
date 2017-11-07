%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function W = wavelet_basis(xsize)

% Define W*x and W'*y for a wavelet transform W
% Input:
%        xsize = a vector of 2 positive integers
% Output:
%        W = struct of 2 fields
%            1) W.times: W*x
%            2) W.trans: W'*y

W.times = @(x) Wavedb1(x,xsize,0); 
W.trans = @(x) Wavedb1(x,xsize,1); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = Wavedb1(x,xsize,trans)
persistent s;
level = 2;
wname = 'db1';

if ~trans;  % W*X
    [Y,s] = wavedec2(reshape(x,xsize),level,wname);
else        % W'*X
    if ~exist('s','var')
    [~,s] = wavedec2(reshape(x,xsize),level,wname);
    end
    Y = waverec2(x,s,'db1');
end
y = Y(:);
 
