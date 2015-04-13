function p = objective(X_1, g_2, X_2, g_3, X_3)
    p = norm(X_1,'fro').^2 + g_2*norm(X_2(:),1) + g_3*norm(svd(X_3),1);
end