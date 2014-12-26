function U_out = updateSubspace(U_in, residual, w, param)
%  update U
   
    eta = param.eta;    
    residualNorm = norm(residual);
    wNorm  = norm(w);
    
    sigma = param.lambda * residualNorm * wNorm;
    
    p = residual / residualNorm;
    
    q = w / wNorm;

  %  t = eta * sigma; % dynamic size

     t =eta;
    
    U_out = U_in + (cos(t)-1) * U_in * (q * q')  - sin(t) * p * q';
    
%       normU = norm(U_out)
%       diffI = norm(U_out'*U_out - eye(size(U_out,2)))
end


 