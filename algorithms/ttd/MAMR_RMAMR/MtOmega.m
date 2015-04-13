function Y=MtOmega(X,I,J,mm,nn)
% Copyright:Xinchen YE, Tianjin University, 2014
    SIZ=[mm nn];
    tmpX = X(:);
    idx = sub2ind(SIZ,I,J);
    tmpY = tmpX(idx);%%%%%%%%%%%%%%%%%%%%%%%%
    %Y=spconvert([I J tmpY;m n 0]);
    Y = sparse(I,J,tmpY,mm,nn);
