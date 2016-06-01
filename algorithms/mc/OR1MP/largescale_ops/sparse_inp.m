% Inner-product and apply sparse selection                   
% This function performs the following:
%       UV = U * V, out(ii) = UV(iR(ii), iC(ii)); 
% To use this function
%     out = sparse_inp_native (U', V, iR, iC, length(iR))
%  U: m x r (for C simplicity, we follow LMaFit to use U')
%  V: r x n
%  iR: len x 1
%  iC: len x 1
%  out: len x 1