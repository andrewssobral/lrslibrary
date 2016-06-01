function [ out ] = sparse_inp_native( U, V, iR, iC )
%SPARSE_INP inner-product and apply sparse selection.
%   Perform the following:
%       UV = U * V, out(ii) = UV(iR(ii), iC(ii)); 

    len = length(iC);
    out = zeros(len, 1); 
    
    for ii = 1: len
        out(ii) = U(iR(ii), :) * V(:,iC(ii)); 
    end

end

