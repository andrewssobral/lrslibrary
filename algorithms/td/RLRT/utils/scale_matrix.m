function X = scale_matrix( X, s, mode )
% X = rectangular matrix whose rows or cols are to be scaled
% s = diag of the scaling diag matrix
% mode = 0  row scaling
%      = 1  col scaling

% prod = sparse( size(X,1), size(X,2) );
if mode
%     for j = 1:size(X,2)
%         X( :, j ) = X(:,j) * s(j);
%     end
    X = scale_cols( X, s );
else
%     for i = 1:size(X,1)
%         X( i, : ) = X(i,:) * s(i);
%     end
    X = scale_rows( X, s );
end
end